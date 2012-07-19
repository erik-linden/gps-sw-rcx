function [cst_kp, w_if_kp, w_df_err_kp, phi_if_kp, isTracking, chip_rate_kp] = PREDICT(...
    prn, lock, magnitude, cst_k, w_if_k, phi_if_k, cst_err_k, phi_err_k, phi_act_k)
% [cst_kp, w_if_kp, phi_if_kp, isTracking, w_df_err_kp] = PREDICT(...
%    prn, lock, P, cst_k, w_if_k, phi_if_k, cst_err_k, phi_err_k, phi_act_k)

CONSTANTS_H;

persistent phi_err_km phi_err_km2 w_df_uneObs_km magnitudeSmooth
persistent w_df_dot_k phi_act_km trackHist chipRateCorr

% If this is the first time we run, initalize everything to zero.
if INITIALIZE
    phi_err_km  = 0;    %previous phase error
    phi_err_km2 = 0;    %next previous phase error
    w_df_uneObs_km = 0; %previous unexplained doppler
    magnitudeSmooth = magnitude;%smoothed magnitude
    phi_act_km = 0;     %previous phase angle
    w_df_dot_k = 0;     %current frequency drift rate
    trackHist = LOL_M;  %number of valid data points in previous M points
    chipRateCorr = 0;   %chip rate correction, used in DLL
end

if USE_AIDING
    % Get predictiv data relevant for tracking:
    % - Predicted apperent Doppler frequency
    % - Predicted actual Doppler frequency
    % - Predicted apperent change in Doppler frequency
    % - Predicted actual change in Doppler frequency
    % - Predicted mean time derivative of Doppler frequency phase
    [w_df_expObs, ~, delta_w_df_expObs, ~, ~] = GET_TRACKING_DATA(cst_k, prn);
else
    w_df_expObs = 0;
    delta_w_df_expObs = 0;
end

% Observed doppler frequency.
w_df_actObs = w_if_k - W_FC;

% Difference in expected apperent Doppler and observed Doppler.
w_df_uneObs = w_df_actObs - w_df_expObs;            

if lock
    % If curve fitting managed a lock, update noise estimates and
    % calculate SNRs.
    magnitudeSmooth = magnitudeSmooth*(1-MAGNITUDE_SMOOTHING) + magnitude*MAGNITUDE_SMOOTHING;
    
    % Mean SNR.
    snrSmoothDb = SNR(magnitudeSmooth);
    
    % CNo ratio, used to optimize bandwidth.
    cnoRatio = 10^(CNO(magnitudeSmooth)/10);
    
    % Instantaneus SNR.
    snrInstDb = SNR(magnitude);

    if snrSmoothDb > DB_LIMIT_SMT && snrInstDb > DB_LIMIT_INS
        % If both mean SNR and instantaneus SNR are above there respective 
        % limits, the data point is valid.
        trackHist = min(trackHist+1,LOL_M);
        
        if trackHist >= LOL_N
            isTracking = true;
            updateErrors = true;
        else
            isTracking = false;
            updateErrors = true;
        end
    else
        % At least one SNR was to low, so the data point is not valid.
        trackHist = max(trackHist-1,0);

        isTracking = false;
        updateErrors = false;
    end
else
    % The curve fitting didn't lock, so the data point is not valid. 
    % We update the magnitude estimate with noise.
    if ~isnan(magnitudeSmooth)
        magnitudeSmooth = magnitudeSmooth*(1-MAGNITUDE_SMOOTHING);
    else
        magnitudeSmooth = 10^(DB_LIMIT_SMT/20);
    end
    trackHist = max(trackHist-1,0);
    isTracking = false;
    updateErrors = false;
end

if isTracking
    % We can track using the current data point.
        
    % Select tracking method. PLL_3RD and FLL updates internal error
    % histories. This gives a small error since it isen't controlled by
    % updateErrors, but this is negliable.
    if USE_PLL && cst_k > PLL_SWITCH_TIME
        if PLL_LOOP_ORDER == 2
            delta_w_df_estObs = PLL_2ND(cnoRatio, phi_err_k, phi_err_km);
        else
            delta_w_df_estObs = PLL_3RD(cnoRatio, w_df_uneObs, w_df_uneObs_km, phi_err_k, phi_err_km, phi_err_km2);
            w_df_uneObs_km = w_df_uneObs;
        end
    else
        [delta_w_df_estObs, w_df_dot_kp] = FLL(w_df_dot_k, phi_act_k, phi_act_km);
        w_df_dot_k = w_df_dot_kp;
    end
    
    % I use this exponential filter to generate DLL corrections. It is a
    % bit strange, but seems to work well.
    chipRateCorr = (1-HNUM)*chipRateCorr + HNUM*cst_err_k;

    % Using actual doppler might seem right here, but almost all of the LO
    % offset comes from sampling clock offset, which effects the appearent
    % chip rate, so in the end, this turns out to be more accurate.
    chip_rate_kp = CHIP_RATE*(1 + chipRateCorr + ...
       (w_df_actObs + delta_w_df_expObs + delta_w_df_estObs)/(2*pi*L1));

    tau = CHIPS_PER_CODE/chip_rate_kp;

    w_if_kp = w_if_k + delta_w_df_expObs + delta_w_df_estObs;
    
    phi_if_kp = mod(phi_if_k + tau*w_if_kp, 2*pi);
else
    % See comment above.
    chip_rate_kp = CHIP_RATE*(1 + chipRateCorr + ...
        (w_df_actObs + delta_w_df_expObs)/(2*pi*L1));

    tau = CHIPS_PER_CODE/chip_rate_kp;

    if trackHist >= LOL_X
        % At least X of the previous data points were valid, so do a
        % differential update.
        w_if_kp = w_if_k + delta_w_df_expObs;
                
        phi_if_kp = mod(phi_if_k + tau*w_if_kp, 2*pi);
    else
        % Fewer than X of the previous data points were valid, so set the
        % IF to the expected IF, to begin searching.
        w_if_kp = W_FC + w_df_expObs + delta_w_df_expObs;
                
        phi_if_kp = mod(phi_if_k + tau*w_if_kp, 2*pi);
    end
end

% Compute the frequency error, (w_df_kp-expected - w_df_kp-actually_set)
w_df_err_kp = (w_if_kp - W_FC) - (w_df_expObs + delta_w_df_expObs);

% Finally, set the CST...
cst_kp = cst_k + tau;

% Uppdate error histories, if we should.
if updateErrors
    phi_err_km  = phi_err_k;
    phi_err_km2 = phi_err_km;
    phi_act_km  = phi_act_k;
else
    phi_err_km  = 0;
    phi_err_km2 = 0;
    phi_act_km  = 0;
end