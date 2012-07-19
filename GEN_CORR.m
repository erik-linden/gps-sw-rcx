function corr = GEN_CORR(code)
% corr = GEN_CORR(code)
%
% Generates correlators according to settings in CONSTANTS.
% 
% Input:
%     code  - [1x1023] CA code from CACODEGN
%     
% Output:
%     corr  - [COH_INT_SAMPxN_CORR] matrix of upsampled correlators

CONSTANTS_H;

index = -(N_CORR-1)/2:(N_CORR-1)/2;     %index vector

time = MAX_SPACE*index/(CHIP_RATE*(N_CORR-1));  %time offset for each correlator

corr = zeros(N_CORR,COH_INT_SAM);  %initialize

for n = 1:N_CORR    %for each correlator
    corr(n,:) = DIGITIZE_CA(code,-time(n),COH_INT_SAM);  %note minus on time
end