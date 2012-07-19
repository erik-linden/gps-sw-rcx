function snr = SNR(magnitude, varargin)
% snr = SNR(magnitude, varargin)
% 
% Returns estimated SNR. The extra argument is an optial coherent
% integration time. If it is left unspecified, the global COH_INT_TIME is
% used.
%
% For details on the algorithm, see:
% http://books.google.se/books?id=stTSHdFhrFUC

global COH_INT_TIME

% Check if T is given.
if isempty(varargin)
    T = COH_INT_TIME;
else
    T = varargin{1};
end

% Compute SNR.
snr = 20*log10(magnitude/noise_std(T));