function [carrierPhaseRate, refTime, weights] = DOPPLER(sv)
% [carrierPhaseRate, refTime, weights] = DOPPLER(sv)
% 
% Returns the doppler shift (called carrier phase rate, ie rad/sec).
% refTime is when the measurments were taken, and weights are the quality
% of the measurments.

global DOPPLER_DB_LIMIT DOPPLER_MAX_PHI_ERR COH_INT_TIME DOPPLER_EXL_RNG

% Load the first file to see how big it is.
nSv = length(sv);
[~ , ~, cst, ~, ~, ~, ~ , ~ , ~ , ~, ~] =  LOAD_TRACK(sv(1));
dim1 = length(cst);

% Initialize vectors.
carrierPhaseRate = nan(dim1,nSv);
refTime = nan(dim1,nSv);
weights = zeros(dim1,nSv);

for n = 1:nSv
    [~ , ~, cst, ~, magnitude, phi_err, ~ , ~ , ~ , ~, w_df] =  LOAD_TRACK(sv(n));
        
    % Set the reference time to the center of the coherent integration interval.
    refTime(:,n) = cst + COH_INT_TIME/2;
       
    % The carrier phase rate is the NCO, approximatly
    carrierPhaseRate(:,n) = w_df;
    
    % 1/mag.^2 is the variance for one sample, so for two...
    weights(:,n) = [0;(magnitude(1:end-1).^2/COH_INT_TIME)];
    
    % Exlude data with low CNO.
    indCno = SNR(magnitude)<DOPPLER_DB_LIMIT;
    
    % Exlude data with large phase error.
    indPhi = abs(phi_err)>DOPPLER_MAX_PHI_ERR*pi/180;
    
    % Exclude data where tracking did not lock.
    indLock = isnan(cst);

    % Exclude a range around bad data.
    ind = indCno|indPhi|indLock;
    ind = conv(double(ind),ones(DOPPLER_EXL_RNG,1),'same');
    ind = ind>0.5;
    
    % Set the weight of the bad data to zero.
    weights(ind,n) = 0;
end