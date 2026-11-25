%% ========================================================================
%  LEVEL 1  --  channel delays: when each pulse actually arrives
%
%  Every channel -- green laser, MW, digitizer trigger -- has its own propagation delay
%  through cables, switches and amplifiers. calibChannelDelays measures them by fitting
%  ONE function (Heaviside steps with exponential responses) to the whole fluorescence
%  trace of a staircase sequence, so the pulse generator can compensate.
%
%  FIELD OFF throughout. These are timing measurements and the field only adds structure
%  to the spectrum that nothing here looks at.
%
%  THE ORDER OF THE CHANNEL LIST decides which edges can be seen at all. Microwave reaches
%  the NVs only while its switch is open AND the AWG is playing, so:
%     {greenLaser, MW, trigger}   the switch is already open: the dip is the AWG itself,
%                                 and this MEASURES THE TRIGGER (on, and the end of the
%                                 AWG waveform, which must agree with it)
%     {greenLaser, trigger, MW}   the AWG is already playing: the dip is the switch,
%                                 and this MEASURES THE MW SWITCH (on and off)
%  The laser edges appear in every run, which is the cross-check: they must agree.
%
%  channelDelay must be longer than the slowest response: >= 3 us with the trigger, and
%  10 us lets the MW recovery (~0.5 us fast + several us slow) finish inside its step.
%
%  Delays are read at the 90% point of each edge (the convention setupInfo.json uses);
%  pass o.delayAt = 'onset' for the switch times instead.
%
%  Writes nothing. Each call prints measured vs stored and what to change.
%  Refit a run without re-measuring:  calibChannelDelays([], struct('data', r2))
%
%  1.6 (calibPhotodiodePair) is the two photodiodes against each other: their timing, and
%  whether dividing by the laser monitor actually lowers the noise. Laser only.
%  ========================================================================

%% 1.0  shared settings
o = struct('channelDelay', 20, 'offset', 1, 'repeats', 100, 'averages', 100);

h = getObjByName('Helmholtz');
h.control('off');


%% 1.1  laser alone -- the reference
r1 = calibChannelDelays({'greenLaser'}, o);


%% 1.2  laser, MW, trigger  --  measures the TRIGGER
r2 = calibChannelDelays({'greenLaser','MW','trigger'}, o);


%% 1.3  laser, trigger, MW  --  measures the MW switch
r3 = calibChannelDelays({'greenLaser','trigger','MW'}, o);


%% 1.4  laser, trigger, MW2  --  measures the second MW switch
%
% Setup 7 drives two microwave channels, each through its own amplifier and circulator,
% and they have different path lengths. Their relative delay is calibrated here; their
% relative PHASE at the carrier is a level 4 job (calibMwPair). timeDelay quantises to
% whole AWG samples, which is 64.6 degrees of phase at 2870 MHz.
r4 = calibChannelDelays({'greenLaser','trigger','MW2'}, o);


%% 1.5  cross-checks
%
% results.edges has one entry per commanded edge: .channel, .kind, .commanded, .tSwitch,
% .tSwitchSterr, .visible. The measured delay is tSwitch - commanded.
%
% (a) The laser is measured in every run. A laser delay that changes with the channel list
%     is not a delay, it is something the sequence does to the measurement.
% (b) In 1.2 the end of the dip is the end of the AWG waveform, which LoadAWG starts on
%     the trigger and sizes to the commanded trigger-off. So its "delay" must equal the
%     trigger-on delay. A difference means the waveform length is not what was asked for.
runs  = {r1, r2, r3, r4};
names = {'1.1 laser', '1.2 laser,MW,trigger', '1.3 laser,trigger,MW', '1.4 laser,trigger,MW2'};
delayOf = @(r, ch, kind) edgeDelay(r, ch, kind);

fprintf('\n(a) laser delays in every run (us):\n');
fprintf('    %-24s %-20s %-20s\n', 'run', 'on', 'off');
for k = 1:numel(runs)
    [dOn, eOn]   = delayOf(runs{k}, 'greenLaser', 'on');
    [dOff, eOff] = delayOf(runs{k}, 'greenLaser', 'off');
    fprintf('    %-24s %+.4f +/- %.4f    %+.4f +/- %.4f\n', names{k}, dOn, eOn, dOff, eOff);
end

% Compared at the ONSET of both: the AWG end is a slow recovery and is always read at its
% onset, so comparing it with the trigger-on 90% point would differ by the drop time alone.
[dT, eT] = edgeDelay(r2, 'trigger', 'on',  't0');
[dA, eA] = edgeDelay(r2, 'trigger', 'off', 't0');
fprintf('\n(b) onset: trigger on %+.4f +/- %.4f us, AWG end %+.4f +/- %.4f us  ->  differ by %.1f sigma\n', ...
        dT, eT, dA, eA, abs(dT - dA) / hypot(eT, eA));

fprintf('\nMW  (1.3): on %+.4f, off %+.4f us\n', delayOf(r3, 'MW', 'on'), delayOf(r3, 'MW', 'off'));
fprintf('MW2 (1.4): on %+.4f, off %+.4f us\n\n', delayOf(r4, 'MW2', 'on'), delayOf(r4, 'MW2', 'off'));


%% 1.6  the two photodiodes: fluorescence (input_1) against the laser monitor (input_2)
%
% Laser only, no MW, field off. About 5 minutes. Writes nothing.
%
% PART 1 measures when each photodiode sees the laser edges, both rows from the SAME
% records, so the 0-10 ns record-start ambiguity (the digitizer clock is not locked to the
% PulseStreamer) cancels in their difference. It suggests spcm.channel_delay(2) in
% setupInfo.json -- the per-input read delay the digitizer already supports, never set.
%
% PART 2 decides whether any of it matters, at the readout window the experiments use
% (json detectionDuration and laserInitializationDuration, fixDelays on): the spread of the
% ratio of two identical halves of a repeat -- what a balanced experiment's contrast is built
% from -- with no referencing, with 'referenced' as it is done today (raw division), and with
% both darks subtracted first, each with the monitor window slid +/-5 samples in software.
%
% Read the WHAT IT SAYS block. Its three possible actions, none taken here:
%   channel_delay(2)   hand edit in setupInfo.json, ONLY if part 2 shows the shift pays
%   dark subtraction   a small change in the photodiode class's 'referenced' path, if it pays
%   'nomral'           if referencing raises the noise (the monitor is noisier than the laser)
%
% Re-analyse without hardware:  pd = calibPhotodiodePair(struct('replay', pd));
pd = calibPhotodiodePair();
P1 = setup7Paths();
outDir = fullfile(P1.folder, 'level 1');
if ~exist(outDir, 'dir'); mkdir(outDir); end
save(fullfile(outDir, 'photodiodePair.mat'), 'pd');
fprintf('saved to %s\n', fullfile(outDir, 'photodiodePair.mat'));


function [d, e] = edgeDelay(r, ch, kind, field)
% Delay of one edge: tSwitch (the convention the run used) by default, or 't0' (onset).
    if nargin < 4; field = 'tSwitch'; end
    d = NaN;  e = NaN;
    m = strcmp({r.edges.channel}, ch) & strcmp({r.edges.kind}, kind) & [r.edges.visible];
    if any(m)
        E = r.edges(find(m, 1));
        d = E.(field) - E.commanded;
        e = E.([field 'Sterr']);
    end
end
