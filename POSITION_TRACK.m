function [posHist, gpsTimeHist, recTime, dopHist, resHist] = POSITION_TRACK(sv, r, useIcp)
% [posHist, gpsTimeHist, recTime, dopHist, resHist] = POSITION_TRACK(sv, r)

CONSTANTS_H
CONSTANTS

% Get all pseudoranges.
[rawPseudoRange refTime prRefTime rawWeights] = PSEUDO_RANGE(sv);

if useIcp %use ICP
    
    % Get ICP pseudoranges.
    [icpPseudoRange, icpRefTime] = ICP_PSEUDO_RANGE(sv);
        
    % Get the size.
    nSv = size(icpPseudoRange,2);
    dim1 = size(icpPseudoRange,1);
    dim2 = floor(dim1/r);
    
    icpPseudoRangeResmpld = nan(size(icpPseudoRange));
    for n = 1:nSv
        indGood = ~isnan(icpPseudoRange(:,n));
        indGoodRaw = ~isnan(rawPseudoRange(:,n));
        icpPseudoRangeResmpld(indGoodRaw,n) = ...
            interp1(icpRefTime(indGood,n),icpPseudoRange(indGood,n),...
            refTime(indGoodRaw,n),...
            'linear','extrap');
    end
    
    % Find the ICP bias, using weighting.
    icpBias = nansum((-rawPseudoRange+icpPseudoRangeResmpld).*rawWeights);
    icpBias = icpBias./nansum(rawWeights.*~isnan(rawPseudoRange+icpPseudoRangeResmpld));
    
    % Resample the ICP ranges onto a uniform vector. The first sample might
    % not be from second 0.
    unRange = zeros(dim2,nSv);
    tOffset = floor(min(nanmin(refTime)));
    recTime = r*((1:dim2)'-1/2)*1e-3 + tOffset;

    for n = 1:nSv
        indGood = ~isnan(icpPseudoRange(:,n));
        unRange(:,n) = interp1(icpRefTime(indGood,n), icpPseudoRange(indGood,n), recTime, 'linear', 'extrap');
        unRange(:,n) = unRange(:,n) - icpBias(n); 
    end
    
    % Time-uniform weighting is the most logical.
    smoothWeights = repmat(nansum(rawWeights),dim2,1);
    
else %don't use ICP.
    
    % Resample onto a uniform vector.
    [unRange, smoothWeights, recTime] = RESAMPLE_PSEUDORANGE(rawPseudoRange, refTime, rawWeights, r);
end

% Check how long our data set is.
dim1 = length(recTime);

% An inital position guess.
guess = 1.0e+006*[3.0989 1.0110 5.4639];

% Initialize vectors.
posHist = nan(dim1,3);
gpsTimeHist = nan(dim1,1);
resHist = nan(dim1,length(sv));
dopHist = nan(dim1,4);

% Tell user if we are not using IONEX data.
load(RECEIVER_FILE);
if Rec.IonType==0
    useIon = 0;
    fprintf('Info: Ion corrections are OFF.\n')
else
    useIon = 1;
    fprintf('Info: Ion corrections are ON.\n')
end

% This can take some time, so a waitbar is nice.
tStart = tic;
h = waitbar(0,sprintf('Creating navigational solution...\n'));
for n = 1:dim1
    indGood = ~isnan(unRange(n,:)) & (smoothWeights(n,:)~=0);
    
    if sum(indGood)>=4
        % Pick out the data for the current point on the grid.
        svTimeAtRef = prRefTime - unRange(n,indGood) + recTime(n);
        prAtRef = unRange(n,indGood);
        weightsAtRef = smoothWeights(n,indGood);
        svAvail = sv(indGood);
        
        % Solve for the position.
        [posObs, gpsTimeOffset, dop, res] = SOLVE_POS(svAvail,prAtRef,svTimeAtRef,weightsAtRef,guess,useIon);
        
        % Use the latest solution as the next guess.
        guess = posObs;
        
        % Update history vectors.
        posHist(n,:)       = posObs;
        gpsTimeHist(n)     = prRefTime + recTime(n) + gpsTimeOffset;
        dopHist(n,:)       = dop;
        resHist(n,indGood) = res;
        
        % Don't reload eph. data.
        INITIALIZE = 0; %#ok<NASGU>
    end
    
    % Update the waitbar every 50 sample.
    if ~mod(n,50)
        tElapsed = toc(tStart);
        waitbar(n/dim1,h,sprintf('Creating navigational solution...\nTime left: %.0f sec',tElapsed*(dim1/n-1)))
    end
end

% Clean up the waitbar.
delete(h)