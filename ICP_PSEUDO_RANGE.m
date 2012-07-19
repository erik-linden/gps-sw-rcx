function [icpPseudoRange, refTime] = ICP_PSEUDO_RANGE(sv)
% [icpPseudoRange, refTime, carrierPhaseRate] = ICP_PSEUDO_RANGE(sv)
% 
% Returns integrated carrier phase pseudorange, at the times speciefied in
% refTime.

global L1

nSv = length(sv);
load(sprintf('tracking_hist_%d',sv(1)))
dim1 = length(cst_hist);

icpPseudoRange = nan(dim1,nSv);
refTime = nan(dim1,nSv);

for n = 1:nSv
    load(sprintf('tracking_hist_%d',sv(n)))
    
    % Clean up data vectors.
    indGood = ~isnan(cst_hist);
    cst_clean = cst_hist(indGood);
    w_df_clean = w_df_hist(indGood);
    
    % Define some stuff...
    deltaT = diff(cst_clean);
        
    % Read out the NCO rate.
    ncoRate = w_df_clean(1:end-1);
    
    % Integrate the carrier phase rate.
    icp = [0;cumtrapz(ncoRate.*deltaT)];
    
    % Get time vector.
    icp_t = cst_clean;
    
    % Convert to distance and add to output.
    icpPseudoRange(indGood,n) = -icp/(2*pi)/L1;
    icpPseudoRange([true; ~indGood(1:end-1)],n) = nan;
    refTime(indGood,n) = icp_t;
    refTime([true; ~indGood(1:end-1)],n) = nan;
end