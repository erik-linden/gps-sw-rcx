function [smoothRange, smoothWeights, recTime] = RESAMPLE_PSEUDORANGE(rawRange, refTime, rawWeights, r)
% [smoothRange, smoothWeights, recTime] = RESAMPLE_PSEUDORANGE(rawRange,
% refTime, rawWeights, r)

% Prepair data vectors.
nSv = size(rawRange,2);
dim1 = length(rawRange);
dim2 = floor(dim1/r);
smoothRange = nan(dim2,nSv);
smoothWeights = nan(dim2,nSv);

% The first sample might not be from second 0.
tOffset = floor(min(nanmin(refTime)));
recTime = r*((1:dim2)'-1/2)*1e-3 + tOffset;

% Loop through SVs.
for n = 1:nSv
    % Loop through segments.
    for k = 1:dim2
        % Decide where we expect to find the datapoints, and search 49
        % samples to either side.
        lowerBound = max(r*(k-1)-49,0);
        upperBound = min(r*k+49,dim1-1);
        searchRange = 1+(lowerBound:upperBound);
        
        % Find the first and last samples with CST in the current interval.
        indFirst = find(refTime(searchRange,n) >= r*(k-1)*1e-3 + tOffset,1,'first');
        indLast  = find(refTime(searchRange,n) <= r*k*1e-3 + tOffset,1,'last');
        
        % Adjust for the first index.
        indFirst = indFirst + lowerBound;
        indLast = indLast + lowerBound;
        
        % Expand the selection two adjacent samples.
        indFirst = max(indFirst-1,1);
        indLast = max(indLast+1,1);
        indVec = indFirst:indLast;
        
        % Check how many are valid.
        indGood = ~isnan(rawRange(indVec,n))&(rawWeights(indVec,n)~=0);
        indVec = indVec(indGood);
        nGood = length(indVec);
        
        % If there are less than two, we can't do curve fitting.
        if nGood<2
            yStar = NaN;
            weight = 0;
        else
            % Pick out the current data points.
            y = rawRange(indVec,n);
            x = refTime(indVec,n);
            w = rawWeights(indVec,n);
            
            % Shift and scale the data, for better precision.
            mu = sum(x)/nGood;
            sigma = sqrt(sum((x-mu).^2)/nGood);
            x = (x-mu)/sigma;
            
            % Solve the least squares problem. Switch two a cubic polynomial
            % when the number of points are large, to capture second order
            % effects.
            if nGood>20
                X = [ones(nGood,1), x, x.^2];
            else
                X = [ones(nGood,1), x];
            end
            W = diag(w);
            beta = (X'*W*X)\X'*W*y;
            
            % Get the y-value at the desired time,
            xStar = (recTime(k)-mu)/sigma;
            if nGood>20
                yStar = beta'*[1; xStar; xStar^2];
            else
                yStar = beta'*[1; xStar];
            end
            
            % Variance of the mean is 1/r^2*sum(var_i,i=1->r). The  weights
            % are the reprocital of the variance.
            weight = nGood^2./nansum(1./rawWeights(indVec,n));
        end
        
        % Store the results.
        smoothRange(k,n) = yStar;
        smoothWeights(k,n) = weight;
    end
end