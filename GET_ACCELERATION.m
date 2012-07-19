function [a correction] = GET_ACCELERATION(recTime)
CONSTANTS;

% Add Joakim's data ploter function to the path.
addpath(DATA_PLOTTER)
dataPath = MAIN_FILE_NAME;

% Load FFU data set.
ffuData = FFUDataSet;
ffuData = ffuData.load_decoded_ex(dataPath,'B',1);

% Read out accelerometer, gyro and gyroX4.
[data time info]      = ffuData.getADCChannel(1,[2 1 15]);
[gyro timeG infoG]    = ffuData.getADCChannel(1,[13 11 9]);
[gyro4 timeG4 infoG4] = ffuData.getADCChannel(1,[8 7 5]);

% Use high resolution gyre when possible.
indLarge = abs(gyro4)>1;
gyroAdj = gyro4;
gyroAdj(indLarge) = gyro(indLarge);
gyroAdj = gyroAdj*2*pi;

% [rp rpTime ~] = ffuData.getStatusChannel(0,11);

% indRp = find(rp == 1, 1);
% time = time(1,:)-time(indRp);

% Adjust time to start of drop.
t0 = time(1);
time = time - t0;
timeG = timeG - t0;

% Accelerometers are given in g's.
data = data*9.82;

% Figure out how much to decimate the data
r = ceil(nanmean(diff(recTime))/nanmean(diff(time)));

timeI = recTime(2:end-1);

% Resample acceleration.
ax = interp1(time,smooth(data(1,:),r),timeI,'linear','extrap');
ay = interp1(time,smooth(data(2,:),r),timeI,'linear','extrap');
az = interp1(time,smooth(data(3,:),r),timeI,'linear','extrap');

a = [ax, ay, az];

% Compute accelerometer correction.
rVec = R_INS_APC;
rVec = repmat(rVec,1,length(gyroAdj));
corr = [cross(diff(gyroAdj,1,2)./repmat(diff(timeG),[3,1]),(rVec(:,1:end-1)+rVec(:,2:end))/2) zeros(3,1)]...
    +cross(gyroAdj,cross(gyroAdj,rVec));

% Resample correction.
cx = interp1(timeG,smooth(corr(1,:),r),timeI,'linear','extrap');
cy = interp1(timeG,smooth(corr(2,:),r),timeI,'linear','extrap');
cz = interp1(timeG,smooth(corr(3,:),r),timeI,'linear','extrap');

correction = [cx, cy, cz];
