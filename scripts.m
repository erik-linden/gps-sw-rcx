%% Local plot
% dub 3098830.3422 1010993.2627 5463999.7967 ITRF2008
% dub = [3098830.3422 1010993.2627 5463999.7967];
CONSTANTS

load(sprintf('%sposLS.mat',TRACK_DIRECTORY))
w = 1./sum(dopHist(:,1:3).^2,2);
w(~isfinite(w)) = 0;

pos = nansum(posHist.*repmat(w,1,3))./nansum(w);
% pos = nanmean(posHist);

lla = ecef2lla(pos);
fprintf('Mean position:\n')
fprintf('%.6f%c, %.6f%c, %.1f m\n',lla(1),char(176),lla(2),char(176),lla(3))

posHistProj = project_to_surface(posHist,pos);
velHistProj = project_vel(velHist,pos);

figure % Position plot
plot(recTime,posHistProj)
grid on
legend('E/W','N/S','Alt.')
title('GPS position - timeseries')
xlabel('Time [sec]')
ylabel('Relative position [meter]')

figure % Velocity plot
plot(velRecTime,velHistProj)
grid on
legend('E/W','N/S','Vert.')
title('GPS velcoity - timeseries')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

% addpath('C:\Program Files\MATLAB\R2011a\toolbox\googleearth')

%% Top-down plot

uvw = posHistProj;
r = 1;  %downsampling factor

% ind = (uvw(:,1)<2.9)&(uvw(:,1)>2.5)&(uvw(:,2)<0)&(uvw(:,2)>-.3);
% [C,I]=max(uvw(:,3));
% ind(I)=true;
% [C,I]=min(uvw(:,3));
% ind(I)=true;
ind = 1:r:length(uvw);

scatter(uvw(ind,1),uvw(ind,2),[],uvw(ind,3),'.')
grid on
axis equal
c = colorbar;
ylabel(c,'Altitude [meter]')

% title('GPS position - local plane')
xlabel('Relative position E/W [meter]')
ylabel('Relative position N/S [meter]')

% kmlStr = ge_point(lla(1:r:end,2),lla(1:r:end,1),lla(1:r:end,3),...
%   'altitudeMode','clampToGround',...
%   'iconURL','http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png');
% kmlStr = ge_plot(lla(:,2),lla(:,1));
% kmlStr = ge_plot3(lla(1:r:end,2),lla(1:r:end,1),lla(1:r:end,3),'altitudeMode','relativeToGround');
% ge_output('ge_output.kml',kmlStr);

%% Velocity magnitude
r = 10;  %downsampling factor
s = 10; %smoothing factor
ind = 1:r:length(velRecTime);

plot(velRecTime(ind),smooth(sqrt(sum(velHistProj(ind,1:3).^2,2)),s))
axis tight
title(sprintf('Velocity smoothed %u samples',s))
xlabel('Time [sec]')
ylabel('Velocity [m/s]')
xlim([0,max(velRecTime)])

%% LO
r = 10; %downsampling factor
s = 10; %smoothing factor
ind = 1:r:length(velRecTime);

plot(velRecTime(ind),smooth(loHist(ind),s)/(2*pi))
xlim([0,max(velRecTime)])
title(sprintf('LO Frequency smoothed %u samples',s))
xlabel('Time [sec]')
ylabel('LO error [Hz]')

%% Clock error
absTime = gpsTimeHist - min(gpsTimeHist) + min(recTime);
timeErr = absTime - recTime;

p = polyfit(absTime,timeErr,2);

plot(absTime,(timeErr-polyval(p,absTime))*1e9)
title(sprintf('Third-order clock error\nDrift is %.2g ns/s and %.2g ns/s^2',p(2)*1e9,p(1)*1e9))
ylabel('Residual error [ns]')
xlabel('"True" time [sec]')

%% Range residuals
s = 1000; %smoothing factor

plot(recTime,reshape(smooth(resHist,s),size(resHist)))
title(sprintf('Range residuals smoothed %u samples',s))
xlabel('Time [sec]')
ylabel('Range residual [meter]')
legend(num2str(sv'))

%% Velocity residuals
s = 1000; %smoothing factor

plot(velRecTime,reshape(smooth(velResHist,s),size(velResHist)))
title(sprintf('Velocity residuals smoothed %u samples',s))
xlabel('Time [sec]')
ylabel('Range residual [meter]')
legend(num2str(sv'))

%% Wind-up rate
plot(velRecTime,windUp*L1/C)
title('Carrier wind-up')
xlabel('Time [sec]')
ylabel('Wind-up rate [Hz]');

%% Drag coefficent
deltaT = nanmean(diff(recTime));

[gx,gy,gz] = xyz2grav(posHist(2:end-1,1),...
            posHist(2:end-1,2),...
            posHist(2:end-1,3),...
            deltaT);
        
acc = diff(posHist,2)/deltaT^2;
accIner = acc - [gx,gy,gz];

accTime = (velRecTime(1:end-1)+velRecTime(2:end))/2;

% refAre = 0.45;    % Drag stuff
refAre = pi*(58e-3)^2;
dragConst = 2*1.06/(1.225*refAre);
velAdj = (velHist(1:end-1,:)+velHist(2:end,:))/2;
velMag = sqrt(sum(velAdj.^2,2));
accProj = sum(accIner.*velAdj,2)./velMag;
plot(accTime,-dragConst*accProj./velMag.^2)
xlabel('Time [sec]')
ylabel('Drag coefficient');
title('Drag')