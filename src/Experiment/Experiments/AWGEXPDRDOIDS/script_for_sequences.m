% Load the .mat file
data = load('axes_and_angles.mat');
data2 = load('segmentDurs.mat');

% Access the axes and angles
axes = data.axes;
angles = data.angles;
segmentDurs = data2.segmentDurs;



% Now you can use axes and angles in your MATLAB code
