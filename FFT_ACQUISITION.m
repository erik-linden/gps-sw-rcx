function  [df, cst, magnitude] = FFT_ACQUISITION(in_sig, code, F0)
% [df, cst, magnitude] = FFT_ACQUISITION(in_sig, code, F0)
%
% Performes FFT-based acquisition.
% 
% Input:
%     in_sig    - vector to search for code
%     code      - digitized code replica, arbitrary length
%     
% Output:
%     df        - estimated doppler frequency [Hz]
%     cst       - estimated code start time
%     magnitude - largest magnitude detected
    

CONSTANTS_H;

% Initialize vectors.
N = length(code);

nCodes = ceil(N/FSAMP_MSEC);

if nCodes > 40 % = 2 chips
	M = round(19*20*FSAMP_MSEC);
elseif nCodes > 20 % = 1 chips
	M = round(8*20*FSAMP_MSEC);
else	
	M = round(3*20*FSAMP_MSEC);
end

LEN = 2^nextpow2(N+M-1);
y = [code;zeros(LEN-length(code),1)];
Y = fft(y);

% dat = [];      %for debugging
% xs=[];
max_val = [0 0 0];
time = (0:1:M-1)'*TP; % Time vector.


% Search Doppler frequency space
for fd=-FD_SIZE+F0:FREQ_STEP:FD_SIZE+F0

    freq_argument = 2*pi*(FC+fd)*time;  % Frequency argument for down modulation.
    x = in_sig(1:M).*exp(-1i*freq_argument);  % Mix to baseband.
    x = [x;zeros(LEN-length(x),1)];
%     xs = [xs (abs(sum((x))))];
    zt = ifft(fft(x).*conj(Y));
    z = abs(zt(N:M-1)); % Pick out the relevant segment.
    
%     dat = [dat, z];    % for debuging
    [y,i] = max(z);    % find maximum value and index to max value
    
%     i = (z(i-1)*(i-1)+z(i+1)*(i+1))/(z(i-1)+z(i+1));

    %     Find maximum value over frequency search
    if(max_val(2)<y)
        %if so, replace the value
        max_val = [fd, y, i];    %these are the doppler, max magnitude, index
    end
end

% We are looking at the absolute value, but actually we want the in-phase
% component only. This gives a reasonable adjustment to the magnitude.
df = max_val(1);
magnitude = sqrt(max_val(2).^2-noise_std(N/FS)^2);

% The peak might be more than one msec from the start of the file, due to
% noise or bit transitions. But its better to start all tracks from the
% same msec, so mod here.
cst = mod((max_val(3)-1)*TP,T);  % -1 because first sample is at zero lag