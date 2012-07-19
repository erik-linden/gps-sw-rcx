function cno = CNO(magnitude, varargin)
% cno = CNO(magnitude, varargin)
% 
% Returns estimated CNR based on SNR, integration time, sampling rate, IF
% filter bandwidth and correlation losses. See
% http://books.google.se/books?id=stTSHdFhrFUC for bakground.

global COH_INT_TIME CORR_LOSS IF_BANDWIDTH FS

if isempty(varargin)
    T = COH_INT_TIME;
else
    T = varargin{1};
end

cno = SNR(magnitude,T) - (10*log10(FS*T)+CORR_LOSS) + 10*log10(IF_BANDWIDTH);