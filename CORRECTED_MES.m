function [rangeCorr, wRange, recTime,...
    prRefTime, dopplerCorr, wDoppler, svPos, svVel] = CORRECTED_MES(sv,r)
% sv = [2 5 6 10 12 13 16 21 23 25 29 30 31];
% sv = [1 3 6 11 15 18 19 22 28];
% sv = [2 5 7 8 10 13 23];
% r = 100;

CONSTANTS_H
CONSTANTS

load(RECEIVER_FILE);
if Rec.IonType==0
    useIon = 0;
    fprintf('Info: Ion corrections are OFF.\n')
else
    useIon = 1;
    fprintf('Info: Ion corrections are ON.\n')
end

% Get number of satellites.
nSv = length(sv);

% Get all pseudoranges.
[rawRange refTime prRefTime rawWeights] = PSEUDO_RANGE(sv);

% Resample onto a uniform vector.
[range, wRange, recTime] = RESAMPLE_PSEUDORANGE(rawRange, refTime, rawWeights, r);

% Get the doppler.
[rawDoppler, refTime, rawWeights] = DOPPLER(sv);

% Resample onto a uniform vector.
[doppler, wDoppler] = RESAMPLE_DOPPLER(rawDoppler, refTime, rawWeights, recTime);

% Compute the SV times, at transmittion.
svTime = prRefTime - range + repmat(recTime,1,nSv);

% Initialize vectors.
gpsTime = nan(size(svTime));
rangeCorr = nan(size(range));
dopplerCorr = nan(size(doppler));

% N is the number of samples in time.
N = length(recTime);

% Vectors used to store satellite ECEF coordinates and velocities, at the
% time of transmittion.
svPos = nan(N,nSv,3);
svVel = nan(N,nSv,3);

% Load the Antenna Phase Center corrections.
apc = APC_CORR(sv);

% Adjust the raw ranges and doppler measurements.
for n = 1:nSv
    
    % Check if we have any ion data.
    if useIon
        % Get the ion delay and differential-code-bias.
        meanSvTime = nanmean(svTime(:,n));
        guess = GET_REC_DATA(nanmean(recTime),'pos');
        pos = GET_EPH(meanSvTime,sv(n),'pos');
        [t1 t2] = ION_CORR(meanSvTime,guess,pos,sv(n));
        tIon = t1;
        tDcb = t2;
    else
        tIon = 0;
        tDcb = 0;
    end
    
    % Use 0 as an intial guess for the clock correction.
    dClk = 0;
    for i = 1:N
        
        % If the SV time is unknown, there is nothing to correct.
        if ~isnan(svTime(i,n))
            
            % The clock error correction is given in gps-time, as
            %  T_SV - dClk(T_GPS) = T_GPS
            %
            % This is solved by fixed point iteration.
            for k = 1:3
                dClk = GET_EPH(svTime(i,n) - dClk,sv(n),'dClk');
            end
            
            % Compute the GPS-time for the current SV, at the time of transmition.
            gpsTime(i,n) = svTime(i,n) - dClk;
            
            % Compute the position and velocity of the current SV at the time of
            % transmition.
            pos = GET_EPH(gpsTime(i,n),sv(n),'pos');
            vel = GET_EPH(gpsTime(i,n),sv(n),'vel');
            
            % Store position and velocity for future use.
            svPos(i,n,:) = pos;
            svVel(i,n,:) = vel;
            
            % Compute the relativistic correction. TODO: CHECK SIGN
            tRel = -2*pos*vel'/C^2;
            
            % Compute the total time correction.
            tCorr = dClk + tRel + tIon + tDcb;
            
            % Compute the corrected range.
            rangeCorr(i,n) = range(i,n) + tCorr + apc(n)/C;
            
            % There is one fewer doppler sample than range, since the
            % doppler is the mean between two range samples.
            if i ~=N
                
                % Clock carrier correction.
                carCorr = GET_EPH(gpsTime(i,n),sv(n),'dClkdt');
                
                % Corrected doppler
                dopplerCorr(i,n) = (doppler(i,n)/(2*pi*L1) - carCorr)*C;
                
            end
            
            INITIALIZE = false;
        end
    end
end

svVel = (svVel(2:end,:,:)+svVel(1:end-1,:,:))/2;