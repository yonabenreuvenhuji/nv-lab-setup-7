b = AWGExpWeakMeasCorr;
path = 'G:\My Drive\NV Lab\Setup 7\calibrations\calibration results\';
folder = '20250505 - with switches';
currentPath = [path, folder, '\'];
jsonLocation = [currentPath, 'expProperties_WEAK.json'];
b.loadJsonParams('all', jsonLocation);

b.repeats = 1000;
b.averages = 1;
b.amplitude = -10;

b.frequency = 2804;

b.tau = round(0.05 : 0.0125 : 0.5, 3);
b.piTime = 0.2;
%b.inputConfig = 'ref';
b.balancedSequence = 1;

b.run;