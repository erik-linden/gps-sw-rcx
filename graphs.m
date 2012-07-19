%% INITIALIZE
clear globals;
CONSTANTS_H;
CONSTANTS;
sv = [5 7 8 9 11 15 17 18 19 26 27 28];%[3 5 7 8 15 18 19 21 26 27 28];%[2 5 6 10 12 13 16 21 23 25 29 30 31];%[2 5 10 12 21 23 25 29 30 31];[5 7 8 15 19 26 28];
nSv = length(sv);
l = 296999;
color = varycolor(nSv);
c = 299792458;


%% Plot lock and tracking indicators
% Plot when lock failed (lower dot) and/or tracking was suspended (upper
% dot). The legend shows the percentage locked/tracked.
leg = cell(nSv,1);
t = linspace(0,l*1e-3,l);
figure
hold on
for n = 1:nSv
    load(sprintf('tracking_hist_%d',sv(n)))
    
    sumNotTracked = sum(~tracking_hist(1:l));
    sumNotLocked = sum(~lock_hist(1:l));
    
    ax = plot(t(~tracking_hist(1:l)),n*ones(1,sumNotTracked)+0.1,'.','Color',color(n,:));
    plot([t(1),t(~lock_hist(1:l))],n*[1,ones(1,sumNotLocked)],'.','Color',color(n,:))
    
    set(get(get(ax,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

    leg(n) = {sprintf('SV%02i - %4.2f%% - %4.2f%%',sv(n),100*sumNotLocked/l,100*sumNotTracked/l)};
    
    drawnow
end
legend(leg)
xlim([t(1) t(end)])
ylim([0,nSv+1])
title('Lock and tracking indicators')
xlabel('Time - [sec]')
hold off
fclose all;

%% Plot de-trended CST
% This is good for getting an idea of the noise in pseudorange. The legend
% shows the standard deviation.
indSv = [1:nSv];
k = 1000;

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,(l-1)*1e-3,l);
h1=figure();
hold on
h2=figure();
hold on
% h3=figure();
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))
    
    indGood = ~isnan(cst_hist(1:l));

    p = polyfit(t(indGood),cst_hist(indGood)',2);
    
    y_corr = c*(cst_hist(1:l)'-polyval(p,t(1:l)));
    
    set(0,'CurrentFigure',h1)
    cstVec = cst_hist(indGood);
    yVec = smooth(y_corr(indGood),k);
    indVec = 1:250:length(cstVec);
    plot(cstVec(indVec),yVec(indVec),'Color',color(n,:))

    set(0,'CurrentFigure',h2)
    [f,xi] = ksdensity(y_corr);
    plot(xi,f,'Color',color(n,:))
     
%     set(0,'CurrentFigure',h3)
%     Hs=spectrum.welch;
%     y_corr(~indGood)=0;
%     Hpsd = psd(Hs,y_corr,'Fs',1e3);
%     semilogy(Hpsd.Frequencies,Hpsd.Data,'Color',color(n,:))
%     hold on
    
    leg(n) = {sprintf('SV%02i - %2.0f meter',sv(indSv(n)),nanstd(y_corr))};
    
    drawnow
end
set(0,'CurrentFigure',h1)
% legend(leg)
% title(sprintf('De-trended CST, smoothed %i msec',k))
xlabel('Time [sec]')
ylabel('Detrended pseudorange [meter]')
hold off

set(0,'CurrentFigure',h2)
legend(leg)
title('De-trended CST variation PDF')
xlabel('Offset - [meter]')
ylabel('Probability')
xlim([-100,100])
hold off

% set(0,'CurrentFigure',h3)
% legend(leg,'Location','SouthWest')
% title('De-trended CST variation spectrum')
% ylabel('Offset - [meter]')
% xlabel('Fraction of CIT')
% axis tight
% hold off

fclose all;

%% Plot CST error
indSv = [1:nSv];

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,l*1e-3,l);
h1 = figure();
hold on
h2 = figure();
hold on
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))    
    
    indGood = ~isnan(cst_err_hist(1:l));
    
    m = nanmean(cst_err_hist(indGood));
    s = nanstd(cst_err_hist(indGood));
    
    [f,xi] = ksdensity(c*cst_err_hist);
    
    set(0,'CurrentFigure',h1)
    cstVec = cst_hist(indGood);
    cstErrVec = smooth(cst_err_hist(indGood),1000)*1023e3;
    indVec = 1:250:length(cstVec);
    plot(cstVec(indVec),cstErrVec(indVec),'Color',color(n,:))
    
    set(0,'CurrentFigure',h2)
    plot(xi,f,'Color',color(n,:))
    
    %plot(t(indGood),c*cst_err_hist(indGood),sprintf('.%s',color(n)));
        
    leg(n) = {sprintf('SV%02i %4.0f (%.0f)',sv(indSv(n)),c*m*100,c*s)};
    
    drawnow;
end
set(0,'CurrentFigure',h1)
legend(leg)
ylabel('Chips')
xlabel('Time [sec]')
set(0,'CurrentFigure',h2)
% legend(leg)

% title('CST error PDF')
ylabel('Probability')
xlabel('Prediction error [meter]')
xlim([-60,60])
hold off
fclose all;

%% Plot SNR/CNo
% SNR is the power seen by the receiver, after processing gains. The legend
% shows mean SNR/mean CNo. CNo is an estimate of the receive
% Carrier-to-Noise-Density. The actual quota peak/noise is given by
% 10^(SNR/20). SNR > 15 is usually required to succesfully track. CNo is
% typically 40-45.
indSv = [1:nSv];
k = 10;

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,l*1e-3,l);
h1 = figure();
hold on
h2 = figure();
hold on
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))
    t = cst_hist;
    
    indGood = ~isnan(cst_hist(1:l));
 
    mS = mean(SNR(magnitude_hist(indGood)));
    mC = mean(CNO(magnitude_hist(indGood)));
    
    [f,xi] = ksdensity(magnitude_hist(indGood),'support','positive');
%     [f,xi] = hist(magnitude_hist(indGood),100);
    
    set(0,'CurrentFigure',h1)
    cnoVec = CNO(smooth(magnitude_hist(indGood),k));
    tVec = t(indGood);
    indVec = 1:1:length(cnoVec);
    plot(tVec(indVec),cnoVec(indVec),'.','Color',color(n,:));
        
    set(0,'CurrentFigure',h2)
    plot(xi./noise_std(COH_INT_TIME),f,'Color',color(n,:))
    
    leg(n) = {sprintf('SV%02i - %.1f/%.1f',sv(indSv(n)),mS,mC)};
    
    drawnow;
end
set(0,'CurrentFigure',h1)
% legend(leg)
% grid on
% axis tight
ylim([20,50])
% xlim([0,8])
% title(sprintf('SNR/CNo, smoothed %i msec',k))
xlabel('Time [sec]')
ylabel('CN_0')
hold off

set(0,'CurrentFigure',h2)
legend(leg)
title('PDF of signal magnitude')
xlabel('Magnitude - [1/std_{noise}]')
ylabel('Probability')
hold off

fclose all;

%% Plot phase error
% The phase error should not be "wraping around". The legend shows the
% standard deviation.
indSv = [1:nSv];

nInd = length(indSv);
leg = cell(nInd,1);
figure
hold on
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))    
    t = cst_hist;
    
    indGood = ~isnan(phi_err_hist(1:l));
    
    m = std(phi_err_hist(indGood));
    
    plot(t(indGood),180/pi*phi_err_hist(indGood),'.','Color',color(n,:));
        
    leg(n) = {sprintf('SV%02i - %2.0f deg',sv(indSv(n)),m*180/pi)};
end
% legend(leg)
% title('Phase error')
xlabel('Time [sec]')
ylabel('Phase error [degree]')
axis tight
ylim([-90,90])
hold off
fclose all;


%% Plot carrier NCO error
% This is the difference between the frequency the carrier tracking NCO ran
% at and the frequency expeced from the ephemeridis and the receiver data
% files.
indSv = [1:nSv];

nInd = length(indSv);
leg = cell(nInd,1);
figure
hold on
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))
    indGood = ~isnan(cst_hist(1:l));
        
    plot(cst_hist(indGood),w_df_err_hist(indGood)/(2*pi),'.','Color',color(n,:));
        
    leg(n) = {sprintf('SV%02i',sv(indSv(n)))};
end
legend(leg)
title('NCO doppler error')
xlabel('Time - [sec]')
ylabel('Frequency error - [Hz]')
grid on
hold off
fclose all;

%% Plot doppler
% This is the actual carrier doppler, measured by phase observations and
% shifted to zero mean. The smoothing time can be altered.
indSv = [1:nSv];
k = 1;

nInd = length(indSv);
leg = cell(nInd,1);
figure
hold on

[carrierPhaseRate, refTime, weights] = DOPPLER(sv(indSv));
carrierPhaseRate = carrierPhaseRate/(2*pi);
for n = 1:nInd
    indGood = ~isnan(carrierPhaseRate(:,n));
    m = mean(carrierPhaseRate(indGood,n));
    
    plot(refTime(indGood,n),detrend(carrierPhaseRate(indGood,n)),'Color',color(n,:));
    
    leg(n) = {sprintf('SV%02i %+5.0f Hz',sv(indSv(n)),m)};
    
    drawnow;
end
legend(leg)
title(sprintf('Detrended carrier, smoothed %i msec',k))
xlabel('Time - [sec]')
ylabel('Frequency - [Hz]')
hold off
fclose all;

%% Plot integrated carrier phase
indSv = [1:nSv];

nInd = length(indSv);
leg = cell(nInd,1);
figure
hold on

[icpPseudoRange, refTime] = ICP_PSEUDO_RANGE(sv(indSv));
for n = 1:nInd
    indGood = ~isnan(icpPseudoRange(:,n));
    pfit = polyfit(refTime(indGood,n),icpPseudoRange(indGood,n),2);
    icpDet = icpPseudoRange(indGood,n) - polyval(pfit,refTime(indGood,n));
    
    plot(refTime(indGood,n),icpDet*C,'-','Color',color(n,:));
        
    leg(n) = {sprintf('SV%02i',sv(indSv(n)))};
end
legend(leg)
title('Detrended integrated carrier phase')
xlabel('Time - [sec]')
ylabel('Distance - [meter]')
hold off
fclose all;
%% Plot bit polarity
indSv = [1:nSv];

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,l*1e-3,l);
figure
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))
        
    subplot(1,2,1)
    hold on
    plot(t,mod(phi_if_hist(1:l)+pi/2,2*pi),'.','Color',color(n,:));
    
    subplot(1,2,2)
    hold on%,1+0.3*n+randn(l,1)/40
    a = phi_if_hist(1:l);
    m = SNR(magnitude_hist(1:l));
    plot(m.*cos(a),m.*sin(a),'.','Color',color(n,:),'MarkerSize',1)
        
    leg(n) = {sprintf('SV%02i',sv(indSv(n)))};
end
subplot(1,2,1)
legend(leg)
title('Bit polarity')
xlabel('Time - [sec]')
ylabel('Angle - [rad]')
ylim([0,2*pi])
xlim([0,l/1000])
hold off
subplot(1,2,2)
axis square
title('Constallation')
xlabel('I')
ylabel('Q')
ylim([-40,40])
xlim([-40,40])
hold off
fclose all;

%% Plot curve fitting RMSE
indSv = [1:nSv];
k=10;

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,l*1e-3,l);
h1=figure();
hold on
h2=figure();
for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))    
    
    % Use the chi distrubution to estimate the RMS of noise
    
    normRmse = rmse_hist(1:l)./(noise_std(COH_INT_TIME));
    indGood = ~isnan(normRmse);
    
    m = mean(normRmse(indGood));
    
    set(0,'CurrentFigure',h1)
    plot(t(indGood),smooth(normRmse(indGood),k),'Color',color(n,:));
    
    set(0,'CurrentFigure',h2)
    Hs=spectrum.welch;
    normRmse = normRmse - m;
    normRmse(~indGood) = 0;
    Hpsd = psd(Hs,normRmse,'Fs',1e3);
    semilogy(Hpsd.Frequencies,Hpsd.Data,'Color',color(n,:))
    hold on
        
    leg(n) = {sprintf('SV%02i - %2.1f %%',sv(indSv(n)),m*100)};
end

set(0,'CurrentFigure',h1)
legend(leg)
title(sprintf('Curve fitting RMSE/magnitude, smoothed %i msec',k))
xlabel('Time - [sec]')
ylabel('RMSE/magnitude')
hold off

set(0,'CurrentFigure',h2)
legend(leg,'Location','SouthWest')
title('RMSE amplitude spectrum')
ylabel('Error')
xlabel('Frequency - [Hz]')
axis tight
hold off

fclose all;

%% SNR STFT
indSv = [1:nSv];

minF = 0;
maxF = 10;

win = 2^8;%1e3/minF;
deltaF = (maxF-minF)/50;

nInd = length(indSv);
leg = cell(nInd,1);
t = linspace(0,l*1e-3,l);
ha = tight_subplot(1,nInd,0.01,.1,.075);

for n = 1:nInd
    load(sprintf('tracking_hist_%d',sv(indSv(n))))
    
    indGood = ~isnan(cst_hist(1:l));
    
    m = mean(magnitude_hist(indGood));
    magCorr = magnitude_hist-m;
    magCorr(~indGood) = 0;
    
%     subplot(1,nInd,n)
    axes(ha(n));
    spectrogram(magCorr,win,win/2,minF:deltaF:maxF,1E3);
    xlabel('')
    ylabel('')
    title(sprintf('SV%02i',sv(n)))
end
set(ha(2:end),'XTickLabel',''); set(ha(2:end),'YTickLabel','')

xlabel(ha(1),'[Hz]')
ylabel(ha(1),'Time [sec]')
hold off
fclose all;

%% Read out elevation and doppler
elev = zeros(1,nSv);
df = zeros(1,nSv);

for n = 1:nSv
    [dfSv, elevSv] = ESTIMATE_DOPPLER(0, sv(n));
    df(n) = dfSv;
    elev(n) = elevSv;
end

%% Plot 
figure;
nMax = 1;
step = 5;
start = 190;
CONSTANTS

for n = 1:nMax
    subplot(nMax,1,n)
    sec = (n-1)*step+start;
    sig = LOAD_GPS_DATA(FILE_NAME,sec);
    nfft = 2^(nextpow2(length(sig))-1);
    [Pxx,f] = pwelch(sig,nfft,0,FC+(-5e3:100:5e3),FS);
    plot(f,Pxx)
    drawnow
%     ylim([-70,-50])
    title(sprintf('Second %01.0d, mean %.2f %%',sec,100*mean(sig)/std(sig)))
end