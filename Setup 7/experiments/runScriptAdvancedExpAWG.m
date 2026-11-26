..................%% Start working
% path = 'G:\My Drive\NV Lab\Setup 7\old files\old-calibrations\calibration results\';
% folder = '20250505 - with switches';
% currentPath = [path, folder, '\'];
% jsonLocation = [currentPath, 'expProperties.json'];
% h = Helmholtz;
% h.load([currentPath, 'helmholtzState.mat']);
% h.control('auto');
% 
 main; 
P            = setup7Paths; %edited by yarden based on yonas changes to basic scripts
jsonLocation = P.json;
currentPath  = [P.folder filesep];

h = Helmholtz;
h.load(P.json);              % the "Helmholtz" block of expProperties.json (falls back to an old helmholtzState.mat)
h.control('auto');
% 



%% single drive
x0 = 0.016;
fRabi = 10;

s1 = AWGsensingSingleDriveCont;

s1.loadJsonParams('all', jsonLocation);

s1.repeats = 1000;
s1.averages = 1;
s1.amplitude = -22;
s1.rabiFreq1 = fRabi;
s1.signalRatio = -10;

s1.tau = x0 + (0 : 1/fRabi : 1.5);

s1.recordVoltageSpan = 1;

s1.laserInitializationDuration = 15;

%s1.inputConfig = 'norm';
s1.balancedSequence = 1;

s1.run;
[pulses, f, T, coeff, ~, ~] = rabiFit(s1, 0, 0, 1);

% pulses = round(pulses, 3);
% paramsStruct = struct('halfPiTime', pulses(1), 'piTime', pulses(2), 'threeHalvesPiTime', pulses(3));
% JsonInfoReader.setParams(paramsStruct, jsonLocation);

%% Ico type one
w = AWGExpTypeOneBase; %I removed constant time and amplitude as they were defined in the ExperimentAWG yarden 26.8.26
w.isTracking = 0; w.nChannels = 2; w.GIMeas = 0; w.AWGorSRSswitch = 0; w.lastDelay = 1;
w.repeats = 3000;
w.averages = 3;
w.referenceDetectionDuration = 0.5;
w.detectionDuration = 0.5;
w.timeDelay = 20e-6;
w.frequency = 2820.6;%2756;
w.halfPiTime = 0.0539/2;
w.piTime = 0.0539;
normalize = 0.1;    % ?
w.segmentDurs = w.segmentDurs * normalize/0.2; 
%obj.segmentDurs = obj.segmentDurs * normalize/0.2; % obj? changed by
%yarden 26.8.26
w.mormSegmentDuration=normalize;
w.EchoLike = 0;
w.amplitude =-16 + [0 0];% -3 + [0 0];
w.nRepetitions = 1:1:120;
w.RabiFreq = 11.227;    % ?
w.RabiX0 = 0;   % ?

w.run;



%% Ico type droid
dr = AWGExpDroid;
dr.isTracking = 0; dr.nChannels = 2; dr.GIMeas = 0; dr.AWGorSRSswitch = 0; dr.lastDelay = 1;
dr.repeats = 3000;
dr.averages = 3;
dr.referenceDetectionDuration = 0.5;
dr.detectionDuration = 0.5;
dr.timeDelay = 20e-6;
dr.frequency = 2756;
normalize = 10;
dr.segmentDurs = dr.segmentDurs * normalize/0.02;
dr.mormSegmentDuration=normalize;
dr.EchoLike = 0;
dr.amplitude = -3 + [0 0];
dr.nRepetitions = 1:3:50;
dr.RabiFreq = 14.0;
dr.RabiX0 = 0;

dr.run;

%% Multi-pulses dynamical decoupling
p = AWGPulsedDynamicalDecoupling;
p.loadJsonParams('all', jsonLocation);

p.repeats = 500;
p.averages = 10;

p.inputConfig = 'diff';
p.balancedSequence = 1;
p.smallDelay = 0.5;
p.constantTime = 1;

p.tau = 78;

p.cycles = 1:1:10;
p.pulsePhase = 'Y';

p.run;

