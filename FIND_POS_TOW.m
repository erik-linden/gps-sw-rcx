function towVec = FIND_POS_TOW(bit)
% towVec = FIND_POS_TOW(bit)
% 
% Finds tentative time-of-week numbers. The TOW is given relative the
% leding edge of the first bit.

% Correlate with preamble.
c = abs(xcorr(bit(:,1)*2-1,[1 -1 -1 -1 1 -1 1 1]));

% This is the index of the first bit of the preamble. The length of c is
% 2*length(bit)-1).
ind = find(c(1:end-59)==8)-length(bit)+1;

towVec = [];

for n =1:length(ind)
    % Make sure there are no NaNs in the data. The TOW is in the second
    % word and D29 and D30 in the first word are needed for the parity
    % check.    
    anyNans = sum(isnan(bit(ind(n)+(28:59))));
    
    if ~anyNans
		% Pick out the current word, and the two previous bits.
        word = bit(ind(n)+(30:59))';
        D29  = bit(ind(n)+(28));
        D30  = bit(ind(n)+(29));
        
        % This corrects polarity to.
        [data, dataOk] = PARITY_CHECK(word,D29,D30);
        
        % If data passes parity...
        if dataOk
            % TOW is bit 1-17.
            s = int2str(data(1:17));
            
            % Get TOW at first bit. TOW is the transmition of the first
            % bit of the _next_ subframe.
            tow = bin2dec(s)*6-6-(ind(n)-1)*20e-3;
            
            % Add to list.
            towVec = [towVec; tow];
        end
    end
end