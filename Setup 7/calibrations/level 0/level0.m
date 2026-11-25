%% ========================================================================
%  LEVEL 0  --  the optical chain, with no microwave and no field involved
%
%  Two measurements, in this order, because the second needs the first.
%
%    calibDarkLevel        what the photodiode reads with the laser blocked. Every power
%                          number downstream is a difference against this, so an
%                          uncorrected offset propagates into all of them.
%    calibLaserSaturation  fluorescence versus excitation power, and the knee. Sets the
%                          operating power: high enough to polarise the NV in a reasonable
%                          time, low enough that more power buys nothing but heat.
%
%  NOTHING in level 0 touches the AWG or the coils, so it can be run with the microwave
%  chain in any state. Run it first after any optical change -- a moved fibre, a cleaned
%  objective, a re-aligned beam -- because everything above it assumes these two numbers.
%
%  Writes NOTHING to expProperties.json, so this level has no write cell. There is no dark-level
%  or laser-power field there and nothing reads one (both functions say so in their WHAT TO DO
%  blocks). The one actionable number, the AOM setting, lives on the laser object: 0.2 prints
%  the line to set it.
%
%  WHAT GETS PLOTTED. Cell 0.1 draws nothing, deliberately -- it produces three numbers per
%  channel and a plot of three numbers is a worse table. Cell 0.2 opens ONE figure with two
%  panels: fluorescence against measured power with the fit and its linear background, and
%  the AOM voltage against measured power. If the fit could not be trusted it opens a
%  different two-panel figure showing the raw normalised data instead, and says so.
%  ========================================================================

%% 0.1  dark level
%
% Block the laser before running this. It will not check for you, and a dark level measured
% with the laser on is a large positive offset that quietly cancels part of every later
% power reading.
d = calibDarkLevel();

% calibDarkLevel returns the two levels SEPARATELY rather than one number, because the two
% answer different questions and only one of them is the thing you subtract:
%
%   d.digitizerOffset    channel A -- the digitizer alone, photodiode unpowered. This is the
%                        instrument's own zero and is not a property of the optics.
%   d.darkWithDetector   channel B -- detector powered, laser gate closed. THIS is the level
%                        every later fluorescence reading sits on top of, and the one
%                        calibLaserSaturation wants.
%   d.detectorContribution   B - A, from PAIRED cycles rather than the difference of two
%                        separate averages, so slow drift cancels. Consistent with zero
%                        means the photodiode contributes no measurable dark current and no
%                        green is leaking past the gate.
%
% Each of A and B is a per-channel struct (.mean, .sterr, .nRepeats), input_1 being
% fluorescence and input_2 the laser monitor -- hence the indices below.
fprintf('\n--- level 0.1 summary ---\n');
fprintf('dark level to subtract (B, detector powered, gate closed):\n');
fprintf('   fluorescence  %+.6f V  +/- %.6f\n', ...
        d.darkWithDetector.mean(1), d.darkWithDetector.sterr(1));
fprintf('   laser monitor %+.6f V  +/- %.6f   (%d reads over %d cycles)\n', ...
        d.darkWithDetector.mean(2), d.darkWithDetector.sterr(2), ...
        d.darkWithDetector.nRepeats, d.nCycles);
fprintf('detector contribution B - A (should be consistent with zero):\n');
fprintf('   fluorescence  %+.6f V  +/- %.6f\n', ...
        d.detectorContribution(1), d.detectorContributionSterr(1));
fprintf('   laser monitor %+.6f V  +/- %.6f\n', ...
        d.detectorContribution(2), d.detectorContributionSterr(2));

% A dark level is only meaningful if it is small against the signal it will be subtracted
% from. Setup 7's fluorescence runs ~20 mV, so anything past a few hundred microvolts on
% input_1 means the gate is leaking or the laser is not actually blocked.
if abs(d.darkWithDetector.mean(1)) > 5e-3
    warning(['level0: the fluorescence channel reads %.3f mV with the gate closed. That is ' ...
             'a quarter of the working signal. Check the laser is really blocked before ' ...
             'using this as a zero.'], 1e3*d.darkWithDetector.mean(1));
end


%% 0.2  laser saturation
%
% Unblock the laser. The dark level from 0.1 is passed in rather than re-measured, so the
% two share one reference -- measuring it twice would put the difference between two
% measurements into the saturation curve. calibLaserSaturation takes the whole struct from
% 0.1 and picks darkWithDetector out of it itself.
s = calibLaserSaturation([], [], [], d);

% calibLaserSaturation now ends with its own "WHAT TO DO WITH THIS" block, so this cell
% only restates the single number you act on. It deliberately does NOT repeat Psat/Fsat:
% in the linear regime those are degenerate, and printing them again here was how a
% meaningless "Psat 0.025 +/- 1.3" ended up looking like a result.
fprintf('\n--- level 0.2 summary ---\n');
if isfield(s, 'recommendedAom') && isfinite(s.recommendedAom)
    fprintf('SET THE AOM TO %.3f V   ->   laser.aom.value = %.3f;\n', ...
            s.recommendedAom, s.recommendedAom);
    if isfield(s, 'plateauAom') && diff(s.plateauAom) > 0
        fprintf('   (flat from %.3f to %.3f V; the midpoint is the robust pick)\n', ...
                s.plateauAom(1), s.plateauAom(2));
    end
else
    fprintf('no usable brightest setting -- see the messages above.\n');
end
if isfield(s, 'linearRegime') && s.linearRegime
    fprintf('Psat was NOT measurable: no saturation over the available power. Record nothing.\n');
elseif isfield(s, 'Psat') && isfinite(s.Psat)
    fprintf('Psat  %.6f +/- %.6f  (monitor volts)\n', s.Psat, s.PsatSterr);
end
% Drift, quoted against the scan's full swing rather than against the re-measured point --
% a randomly-chosen near-dark first point once turned a 25 uV change into "35.3% drift".
if isfield(s, 'driftSigma')
    fprintf('drift over the scan: %+.1f sigma', s.driftSigma);
    if isfield(s, 'driftPercent')
        fprintf('  (%+.1f%% of the fluorescence swing)', s.driftPercent);
    end
    if isfield(s, 'driftPointUsable') && ~s.driftPointUsable
        fprintf('   [weak check: the re-measured point was near dark]');
    end
    fprintf('\n');
end
