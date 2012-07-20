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
s = 100; %smoothing factor
ind = 1:r:length(velRecTime);

plot(velRecTime(ind),smooth(sqrt(sum(velHistProj(ind,1:3).^2,2)),s))
axis tight
title('Velocity')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')
xlim([0,max(velRecTime)])


%% LO
r = 10;  %downsampling factor
s = 100; %smoothing factor
ind = 1:r:length(velRecTime);

plot(velRecTime(ind),smooth(loHist(ind),s)/(2*pi))
xlim([0,max(velRecTime)])
titel(sprintf('LO Frequency smoothed %u samples',s))
xlabel('Time [sec]')
ylabel('LO error [Hz]')

%% Clock error
absTime = gpsTimeHist - min(gpsTimeHist) + min(recTime);
timeErr = absTime - recTime;

p = polyfit(absTime,timeErr,2);

plot(absTime,(timeErr-polyval(p,absTime))*1e9)
title(sprintf('Third-order clock error\nDrift is %.1G ns/s and %.1G ns/s^2',p(2)*1e9,p(1)*1e9))
ylabel('Residual error [ns]')
xlabel('"True" time [sec]')
%%
k = 1000;
plot(recTime,reshape(smooth(resHist,k),size(resHist)))
xlabel('Time [sec]')
ylabel('Range residual [meter]')
xlim([0,92])

%%
k = 1000;
plot(reshape(smooth(velResHist,k),size(velResHist)))

%%
figure
Hs=spectrum.welch('Hamming',2^11,50);
Hpsd = psd(Hs,posHistProj(:,1),'FreqPoints','User Defined','FrequencyVector',3:.01:6,'Fs',1e3);
% Hpsd2 = psd(Hs,posHistProj(:,2),'FreqPoints','User Defined','FrequencyVector',3:.001:6,'Fs',1e3);
% plot(Hpsd.Frequencies,sqrt(Hpsd.Data.^2+Hpsd2.Data.^2))
plot(Hpsd.Frequencies,Hpsd.Data)


%%
nfft = 2^19;
P = fft(posHistProj(:,1),nfft);
plot(linspace(0,1000,nfft),2*abs(P)/length(posHistProj))


%% Plot velocity vs. time
lla = ecef2lla(velHist);
klm = project_vel(velHist,pos);

plot(velRecTime,klm)
grid on

title('Velocity')
legend('E/W','N/S','Up/Down')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')

%%
indGood = ~isnan(velHist(:,1));
dim1 = length(velRecTime);

velPos = zeros(dim1+1,3);
for n = 1:3
    velClean = interp1(velRecTime(indGood),velHist(indGood,n),...
        velRecTime,'linear','extrap');
    
    tSpl = (velRecTime(1:end-1)+velRecTime(2:end))/2;
    dt = diff(tSpl);
    mDt = mean(dt);
    dt = [dt(1)-mDt, dt', dt(end)+mDt]';
    relPos = [0; cumsum(velClean).*dt];
    
    relPos = relPos - mean(relPos(1:end-10,:)) + nanmean(posHist(1:end-10,n));
    
    velPos(:,n) = relPos;
end

lla = ecef2lla(velPos(1:end-40,:));
% kmlStr = ge_plot(lla(:,2),lla(:,1));
% ge_output('ge_output.kml',kmlStr);

plot(velRecTime,project_to_surface(velPos(1:end-1,:)))

%% 
p = 0.995;

indGood = ~isnan(klm(:,1));

fnplt(fnder(csaps(velRecTime(indGood),klm(indGood,:)',p)))
grid on

%%
ind = ~isnan(velHist(:,1));
temp=fnval(pp,timeVelPos);
col=['b' 'g' 'r'];
for n = 1:3
    plot(recTimeVel(ind),smooth(velHist(ind,n),50))
    hold on
    plot(timeVelPos,temp(n,:),col(n))
end
hold off

%%
% pos = nanmean(velPos);
% 
% uvw = project_vel(velCorr, pos);
% col=['b' 'g' 'r'];
% 
% for n=1:3
%     plot(recTimeVel,smooth(uvw(:,n),50),col(n))
%     hold on
% end
% hold off
% title('Velocity, smoothed')
% legend('E/W','N/S','Vertical')
% xlabel('Time - [s]')
% ylabel('Velocity - [m/s]')
%  
% disp(sqrt(nanmean(uvw(:,1))^2+nanmean(uvw(:,2))^2))
% disp(nanmean(uvw(:,3)))
% 
% %%
% qq = ppdiff(pp,1);
% 
% v=fnval(pp,timeVelPos);
% a=fnval(qq,timeVelPos);
% 
% m = 1;
% A = 1;%pi*0.06^2;
% rho = 1;
% g = 9.82;
% 
% a_in = a;
% a_in(3,:) = a_in(3,:)+g/m;
% 
% cd=-2*dot(v,a_in)*m./(rho*A*sum(v.^2).^(3/2));
% figure
% plot(timeVelPos,cd)
% 
% title('Drag')
% xlabel('Time - [s]')
% ylabel('Drag area per mass - [Cd*A/m]')
% 
% %%
% plot(recTimeVel,smooth(sqrt(sum(velCorr.^2,2)),10))
%     