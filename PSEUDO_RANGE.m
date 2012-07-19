function [pseudoRange, refTime, prRefTime, weight] = PSEUDO_RANGE(sv)
% [pseudoRange refTime prRefTime weight] = PSEUDO_RANGE(sv)

global POSITION_DB_LIMIT POSITION_EXL_RNG COH_INT_TIME

% Load a tracking history file to see how big it is.
load(sprintf('tracking_hist_%d',sv(1)))

nSv = length(sv);
dim1 = length(cst_hist);

% Also find out at what second in receiver time we started. This is
% subtracted from the PR reference time, to keep the PR resonable.
fistSecond = floor(nanmin(cst_hist));

% Compute what the SV clocks were at the transmittion of the first received
% CST.
svTimeAtFirstSam=TOT_FIRST_SAM(sv);

% Compute the correction to the pseudoranges due to differences in SV
% times at the first CSTs.
prRefTime = nanmin(svTimeAtFirstSam) - fistSecond;
totCorrection = svTimeAtFirstSam - prRefTime;

% Initialize vectors.
pseudoRange = zeros(dim1,nSv);
refTime = zeros(dim1,nSv);
weight = zeros(dim1,nSv);

% Compute pseudoranges for each SV.
for n = 1:nSv
    load(sprintf('tracking_hist_%d',sv(n)))
    
    % Doppler shift of the code can cause some cst errors, especially with
    % longer integration times. Here I try to correct for that, by assuming
    % the code will correlate with the middle of the compressed/expanded
    % doppler shifted code. See sketch:
    %
    % signal    |======^=======|  <- %(SV)% (doppler>0)
    % code    |========^=========|
    %         |========| => cst + COH_INT_TIME/2
    
    % Compute the time at the center of the code.
    codeCenterTime = cst_hist + COH_INT_TIME/2;
    
    % Compute pseudorange, relative the "prRefTime" reference time.
    % Start with the code time evolution. If a SV transmitted later than the
    % reference, it is closer, for the same time of reception.
    % Finally, remove the expected increase in range. i.e. 1 msec per code.
    pseudoRange(:,n) = + cst_hist...     
                       - totCorrection(:,n)... 
                       - (0:length(cst_hist)-1)'*1e-3;
                   
    % This is the receiver time when a given range was observed.
    refTime(:,n) = codeCenterTime;
    
    % Signal strength is used to weight the position solution. The weights
    % should be the reprocital of the variance.
    weight(:,n) = magnitude_hist.^2/COH_INT_TIME;
    
    % Exclude data with to low SNR.
    ind = SNR(magnitude_hist)<POSITION_DB_LIMIT;
    
    % Exclude a range around bad data.
    ind = conv(double(ind),ones(POSITION_EXL_RNG,1),'same');
    ind = ind>0.5;
    
    % Too weak signals get weight zero.
    weight(ind) = 0;
end