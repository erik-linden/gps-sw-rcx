function alpha = estimate_noise_alpha(prn, nCodes,file)
% alpha = estimate_noise_alpha(prn,nCodes,file)
%
% Applise the receiver filter to several empty channels, and looks at the 
% variance of the output.

CONSTANTS_H;
CONSTANTS;

in_sig = LOAD_GPS_DATA(RAW_FILE,file);

nSV = length(prn);
alphaN = zeros(1,nSV);
meanVal = zeros(1,nSV);

% Do for all SVs.
for k = 1:nSV
    code = CACODEGN(prn(k));
    
    code=DIGITIZE_CA(code,0,nCodes*ONE_MSEC_SAM);
    
    N = length(code);
    
    LEN = 2^nextpow2(3*N-1);
    x = zeros(3*N-1,1);
    y = zeros(3*N-1,1);
    y(1:N) = code;
    Y = fft(y,LEN);
    
    doppler = linspace(-5e3,5e3,200/nCodes);
    data = zeros(N,length(doppler));      %for debugging
    
    % Search Doppler frequency space
    for n = 1:length(doppler)
        fd = doppler(n);
        
        time = (0:1:2*N-1)'.*TP; % Time vector.
        freq_argument = 2*pi*(FC+fd)*time;  % Frequency argument for down modulation.
        x(N:3*N-1) = in_sig(1:2*N).*exp(-1i*freq_argument);  % Mix to baseband.
        zt = ifft(fft(x,LEN).*conj(Y));
        z = zt(N:2*N-1); % Pick out the relevant segment.
        
        data(:,n) = z;    %for debuging
    end
    
    dataVec = [real(data(:));imag(data(:))];
    
    [muhat,sigmahat] = normfit(dataVec);
    
    a = sigmahat*(nCodes*1e-3)^(-1/2);
    
    alphaN(k) = a;
    meanVal(k) = muhat;
end

alpha = mean(alphaN);

fprintf('Noise parameter alpha is: %.2e\n',alpha)
fprintf('Inter-SV standard deviation is: %.1f %%\n',std(alphaN)/alpha*100)