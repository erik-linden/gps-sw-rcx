function digCode = DIGITIZE_CA(code, dt, N)
% digCode = DIGITIZE_CA(code, dt, N)
%
% Upsamples the CAcode with a time offset. 
% 
% Input:
% 	code	- [1023] CA code from CACODEGN
%  	dt      - time offset at beginning of code
% 	N       - number of samples in the upsampled code
%     
% Output:
% 	digCode - [N] digitized CA code
% 
% dt<0 => sampling starts early
% dt>0 => sampling starts late

CONSTANTS_H;

% Generate a sampling vector.
samp_index = 1:N;

% Generate the time base for sampling from the sampling vector with the
% appropriate dt offset.
time_base = dt + (samp_index-1)*TP;

% Normalize time_bases to range from 0 to n;
time_base = time_base/T;

% Create sample_times which is scaled from 0 to CHIPS_PER_CODE including
% wrap-around.
samp_times = CHIPS_PER_CODE*(time_base-floor(time_base));

% Finally convert to integers for array access.
samp_index = floor(samp_times)+1;

% Replicat the code the specified number of times.
code = repmat(code,ceil(N/ONE_MSEC_SAM),1);

% And up-sample the CACODE.
digCode = code(samp_index,1)*2-1;