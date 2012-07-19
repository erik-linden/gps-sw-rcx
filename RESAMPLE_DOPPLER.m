function [smoothDoppler, smoothWeights] = RESAMPLE_DOPPLER(rawDoppler, refTime, rawWeights, recTime)
% [smoothRange, smoothWeights, recTime] = RESAMPLE_PSEUDORANGE(rawRange,
% refTime, rawWeights, r)

% This is similar to RESAMPLE_PSEUDORANGE, but outputs the mean doppler
% between the times given in the recTime vector. That is also why no cubic
% function is used when there are many data points.

% Prepair data vectors.
nSv = size(rawDoppler,2);
dim1 = length(rawDoppler);
dim2 = length(recTime)-1;
smoothDoppler = nan(dim2,nSv);
smoothWeights = nan(dim2,nSv);

% The first sample might not be from second 0.
tOffset = floor(min(nanmin(refTime)));

% Loop through SVs.
for n = 1:nSv
    % Loop through segments.
    for k = 1:dim2
        % Decide where we expect to find the datapoints, and search 49
        % samples to either side.
        lowerBound = max(floor((recTime(k)-tOffset)/1e-3)-49,0);
        upperBound = min(floor((recTime(k+1)-tOffset)/1e-3),dim1-1);
        searchRange = 1+(lowerBound:upperBound);
        
        % Find the first and last samples with CST in the current interval.
        indFirst = find(refTime(searchRange,n) >= recTime(k),1,'first');
        indLast  = find(refTime(searchRange,n) <= recTime(k+1),1,'last');
        
        % Adjust for the first index.
        indFirst = indFirst + lowerBound;
        indLast = indLast + lowerBound;
        
        % Expand the selection two adjacent samples.
        indFirst = max(indFirst-1,1);
        indLast = max(indLast+1,1);
        indVec = indFirst:indLast;
        
        % Check how many are valid.
        indGood = ~isnan(rawDoppler(indVec,n))&(rawWeights(indVec,n)~=0);
        indVec = indVec(indGood);
        nGood = length(indVec);
        
        % If there are less than two, we can't do curve fitting.
        if nGood<2
            yStar = NaN;
            weight = 0;
        else
            % Pick out the current data points.
            y = rawDoppler(indVec,n);
            x = refTime(indVec,n);
            w = rawWeights(indVec,n);
            
            % Shift and scale the data, for better precision.
            mu = sum(x)/nGood;
            sigma = sqrt(sum((x-mu).^2)/nGood);
            x = (x-mu)/sigma;
            
            % Solve the least squares problem.
            X = [ones(nGood,1), x];
            W = diag(w);
            beta = (X'*W*X)\X'*W*y;
            
            % Get the y-value at the desired time,
            xStar = ((recTime(k)+recTime(k+1))/2-mu)/sigma;
            yStar = beta'*[1; xStar];
            
            % Variance of the mean is 1/r^2*sum(var_i,i=1->r). The  weights
            % are the reprocital of the variance.
            weight = nGood^2./nansum(1./rawWeights(indVec,n));
        end
        
        % Store the results.
        smoothDoppler(k,n) = yStar;
        smoothWeights(k,n) = weight;
    end
end