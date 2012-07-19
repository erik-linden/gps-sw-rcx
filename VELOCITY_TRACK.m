function [velHist loHist recTime dopHist resHist] = VELOCITY_TRACK(sv, posHist, posTime, gpsTime, r)
% [velHist loOffsetHist refTime resHist dopHist] = VELOCITY_TRACK
% (sv,posHist, posTime, gpsTime)

CONSTANTS_H
CONSTANTS

% Get the doppler.
[carrierPhaseRate, refTime, rawWeights] = DOPPLER(sv);

% Resample onto a uniform vector.
[unCpr, weights, recTime] = RESAMPLE_PSEUDORANGE(carrierPhaseRate, refTime, rawWeights, r);

% Get the sizes.
dim1 = length(recTime);
unPos = zeros(dim1,3);

% Resample the position and GPS-time.
indGood = ~isnan(gpsTime);
for n = 1:3
    unPos(:,n) = interp1(posTime(indGood),posHist(indGood,n),recTime,'linear','extrap');
end
unGpsTime = interp1(posTime(indGood),gpsTime(indGood),recTime,'linear','extrap');

% Initialize vectors.
velHist = nan(dim1,3);
loHist = nan(dim1,1);
resHist = nan(dim1,length(sv));
dopHist = nan(dim1,4);

% This can take some time, so a waitbar is nice.
tStart = tic;
h = waitbar(0,sprintf('Creating velocity solution...\n'));
for n = 1:dim1
    indGood = ~isnan(unCpr(n,:));
    
    if sum(indGood)>=4
        % Pick out the data for the current point on the grid.
        gpsAtRef = unGpsTime(n);
        cprAtRef = unCpr(n,indGood);
        posAtRef = unPos(n,:);
        weightsAtRef = weights(n,indGood);
        svAvail = sv(indGood);
        
        % Solve for the velociy.
        [vel, lo, dop, res] = SOLVE_VEL(svAvail, cprAtRef, posAtRef, gpsAtRef, weightsAtRef);
        
        % Update history vectors.
        velHist(n,:) = vel;
        loHist(n) = lo;
        resHist(n,indGood) = res;
        dopHist(n,:) = dop;
        
        % Don't reload eph. data.
        INITIALIZE = 0; %#ok<NASGU>
    end
    
    % Update the waitbar every 50 sample.
    if ~mod(n,50)
        tElapsed = toc(tStart);
        waitbar(n/dim1,h,sprintf('Creating velocity solution...\nTime left: %.0f sec',tElapsed*(dim1/n-1)))
    end
end

% Clean up the waitbar.
delete(h)