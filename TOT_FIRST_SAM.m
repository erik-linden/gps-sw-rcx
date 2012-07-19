function svTimeAtFirstSam=TOT_FIRST_SAM(sv)
% svTimeAtFirstSam=TOT_FIRST_SAM(sv)
%
% Returns the SV time at the transmission of the first CST.

global C

% Read what GPS-week it is.
estimatedStartTime = GET_REC_DATA(0, 'time');

% Initalize vectors.
nSv = length(sv);
svTimeAtFirstBit = zeros(1,nSv);
I = zeros(1,nSv);

% Do for each SV.
for n = 1:nSv
    [~, ~, ~, ~, ~, ~, phi_if , ~, ~, ~, ~] = LOAD_TRACK(sv(n));
    
    % Find the index of the first chip in the first whole bit.
    I(n) = BIT_LOCK(phi_if,sv(n));
    
    % Extract the data bits.
    bit = BIT_EXTRACT(phi_if, I(n));
    
    % Find the TOW at the leading edge of the first whole bit.
    tows = FIND_POS_TOW(bit);
    
    if ~isempty(tows)
        % Convert the TOW to GPS-system time.
        tentativeSvTimes = tows;

        % Find the value that differs the least from the expected value.
        [Y, ind] = min(abs(tentativeSvTimes-estimatedStartTime));
        
        % Check to see if the difference is "to big". Here, 10 minutes.
        if Y < 15*60
            svTimeAtFirstBit(n) = tentativeSvTimes(ind);
        else
            svTimeAtFirstBit(n) = NaN;
        end
    else
        svTimeAtFirstBit(n) = NaN;
    end
end

% Adjust times to get the time at the first CST.
svTimeAtFirstSam = svTimeAtFirstBit-(I-1)*1e-3;

% This is a sort of error checking on the pseudoranges...

% Estimate ranges to SVs.
svRange = ESTIMATE_RANGE(0, sv);

% Compute the difference in expected pseudorange and observed pseudorange.
pseudoRangeDiff = svTimeAtFirstSam+svRange/C;
meanRangeOffset = nanmean(pseudoRangeDiff);
pseudoRangeError = pseudoRangeDiff - meanRangeOffset;

% If the error is more than a chip, throw a warning.
if abs(pseudoRangeError) > 1e-3
    warning('Possible problems when finding Time-of-Week.')
end

% Sometimes there might not be possible to find the TOW for some SVs. In
% that case, it might be possible to guess the correct value.

isGood = ~isnan(svTimeAtFirstBit);
if sum(isGood)<nSv
    fprintf('Info: Some SVs lacks timing, interpolating SV %d.\n',sv(~isGood))
    
    n = round((meanRangeOffset+I*1e-3...
        -svRange/C-max(svTimeAtFirstBit))/20e-3);
    
    svTimeAtFirstBit(~isGood) = n(~isGood)*20e-3+max(svTimeAtFirstBit);
    
    svTimeAtFirstSam = svTimeAtFirstBit-(I-1)*1e-3;
end
