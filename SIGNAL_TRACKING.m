function SIGNAL_TRACKING(doppler_frequency, code_start_time, in_sig, prn, corr, file, fileNo, t0, Nfiles)
% SIGNAL_TRACKING(doppler_frequency, code_start_time, in_sig, prn, corr,
% file, fileNo, t0, Nfiles)
% 
% This is the main tracking loop. Outputs are saved to disk.

CONSTANTS_H;
loadFile = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now initialize variables for parameters and histories       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create storage vectors here, preinitialize for speed

%the following *_hist are perpetual variables that are required to
%determine the navigation solution.  They are stored as
%bit_cst_hist_PRN.mat in the GPS_SW_RCX.m directory

cst_hist        = nan(Nfiles*1000-1,1);   %the code start times
lock_hist       = false(Nfiles*1000-1,1); %lock indicator
magnitude_hist  = nan(Nfiles*1000-1,1);   %signal magnitude
w_df_hist       = nan(Nfiles*1000-1,1);	  %carrier doppler shift (radians)
w_df_err_hist   = nan(Nfiles*1000-1,1);   %difference in expected and set doppler
phi_if_hist     = nan(Nfiles*1000-1,1);	  %carrier phase
tracking_hist   = false(Nfiles*1000-1,1); %tracking indicator
cst_err_hist    = nan(Nfiles*1000-1,1);   %difference in measured cst and NCO cst
phi_err_hist    = nan(Nfiles*1000-1,1);   %difference in measured phi and NCO phi
chip_rate_hist  = nan(Nfiles*1000-1,1);   %NCO chip rate
rmse_hist       = nan(Nfiles*1000-1,1);   %curve fitting RMSE

start = 0; %indexing variable for saving in the above vectors

% Users specify this in relationship to the start of tracking, so we
% simply update the global.
PLL_SWITCH_TIME = PLL_SWITCH_TIME + t0;

% Initialize parameters. Some get values as inputs, others are set to zero.
cst_nco_k       = t0 + code_start_time;
w_df_k          = doppler_frequency*2*pi;
w_df_err_k      = 0;
w_if_k          = W_FC + w_df_k;
phi_if_nco_k    = 0;
isTracking_k    = true;
chip_rate_k     = 1023e3;
index           = 0;

% For longer correlators to work, we need to know when bit transitions are
% going to happen. This is read from the bit-start-time file.
bstExist = exist(BIT_START_TIME_FILE,'file');

if bstExist % The BST file exists, so load it.
    bst_file = load(BIT_START_TIME_FILE,'bst');
    cstLeadingEdgeOfBit = bst_file.bst(prn);
else        % No BST file hase been generated, set all BSTs to zero.
    cstLeadingEdgeOfBit = NaN;
end

% Check if BST data is necissary and if it exists.
if COH_INT_TIME > T && isnan(cstLeadingEdgeOfBit)
    warning('Coherent integration time is >1 msec, but there is no synchronisation data for SV%02d\n', prn)
    
    % Set the current index reletive bit-flips to something arbitrary, like
    % zero.
    indexRelSOB = 0;
else
    % Set the index releative bit-flips. This index is such that index = 0
    % is the first chip in a bit.
    indexRelSOB = mod(round((cst_nco_k-cstLeadingEdgeOfBit)*1000),20);
end

% Stop is the end time of in_sig (minus a buffer).  It tells us we need to
% load the next millisecond of data
stop = floor(length(in_sig)-2*COH_INT_SAM)*TP + t0;

% We now have to iterate over the total number of seconds specified by the
% user.
fileOffset = 0;
tStart = tic;
h = waitbar(0,sprintf('Tracking PRN %02d - XX.X dB, second #%d\n',prn,fileNo));
while(fileOffset < Nfiles)
    if(loadFile)  %if we arrived at the end of the current 1 sec of data
        %load the next second of data
        %update the wait bar
        timeLeft = toc(tStart)*(Nfiles/fileOffset-1);
		timeTotal = toc(tStart)*(Nfiles/fileOffset);
        sigStng = SNR(nanmean(magnitude_hist(start+(1:index))));
        waitbar(fileOffset/Nfiles,h,sprintf('Tracking PRN %02d - %2.1f dB, second #%d\n%d minutes left (%d minutes total)',...
			prn,sigStng,fileNo+fileOffset,ceil(timeLeft/60),ceil(timeTotal/60)));
        loadFile = 0;
        %load the next second of data
        in_sig = LOAD_GPS_DATA(file,fileNo + fileOffset);
        %augment the data with the left over data from the previous second
        in_sig = [samp_buff; in_sig];
        %and update the stop time
        stop = floor(length(in_sig)-2*COH_INT_SAM)*TP + t0;

        %and re-start the index
        start = start+index;
        index = 0;
    end
    
    % Loop through until code_start_time >= stop, which indicates that we
    % have arrived at the end of the current file.
    while(cst_nco_k<stop)

        % Measure current iteration values.
        [lock_k, rmse_k, magnitude_k, cst_k, cst_err_k, phi_if_k, phi_err_k] = MEASURE(...
            in_sig, cst_nco_k, w_if_k, phi_if_nco_k, corr, t0, indexRelSOB);
        
        % Record time k data into history vectors.
        index = index+1;
        cst_hist(start+index)         = cst_k;          %the code start times
        lock_hist(start+index)        = lock_k;         %curve fitting success
        magnitude_hist(start+index)   = magnitude_k;    %signal magnitude
        w_df_hist(start+index)        = w_if_k - W_FC;  %the doppler shift - NCO
        w_df_err_hist(start+index)    = w_df_err_k;     %difference in expected and set doppler
        phi_if_hist(start+index)      = phi_if_k;       %the carrier phase - measured
        tracking_hist(start+index)    = isTracking_k;   %if tracking loop is locking at errors
        cst_err_hist(start+index)     = cst_err_k;      %cst error - difference NCO/measured
        phi_err_hist(start+index)     = phi_err_k;      %phase error - difference NCO/measured
        chip_rate_hist(start+index)   = chip_rate_k;    %chipping rate - NCO
        rmse_hist(start+index)        = rmse_k;         %curve fitting RMSE
        
        % Predict next iteration values.
        [cst_nco_kp, w_if_kp, w_df_err_kp, phi_if_nco_kp, isTracking_kp, chip_rate_kp] = PREDICT(...
            prn, lock_k, magnitude_k, cst_nco_k, w_if_k, phi_if_nco_k, cst_err_k, phi_err_k, phi_if_k);
        
        % Update variables.
        cst_nco_k       = cst_nco_kp;
        w_if_k          = w_if_kp;
        phi_if_nco_k    = phi_if_nco_kp;
        isTracking_k    = isTracking_kp;
        chip_rate_k     = chip_rate_kp;
        w_df_err_k      = w_df_err_kp;
        
        % Increment StartOfBit index.
        indexRelSOB = mod(indexRelSOB+1,20);

        % Unset the initialize flag.
        INITIALIZE = 0;
    end

    % Now save data if code_start_time_k>stop for next iteration of file
    if(cst_nco_k>stop)

        index_left = floor((cst_nco_k-t0)/TP);
        t0 = index_left*TP+t0-TP;

        % Grab what's left of data samples from current time and determine
        % current stop time
        samp_buff = in_sig(index_left:end);
        loadFile = 1;
        fileOffset = fileOffset + 1;
    end
end

% Restore the global to avoid interfering with the next run.
PLL_SWITCH_TIME = PLL_SWITCH_TIME - t0;

% Save all the history vectors.
[status,message,messageid] = mkdir(TRACK_DIRECTORY);
save(sprintf('%stracking_hist_%i',TRACK_DIRECTORY,prn),'cst_hist','lock_hist','rmse_hist',...
    'magnitude_hist','w_df_hist','w_df_err_hist','phi_if_hist','tracking_hist','cst_err_hist','phi_err_hist',...
    'chip_rate_hist','COH_INT_TIME');

% Close the progress bar.
close(h);