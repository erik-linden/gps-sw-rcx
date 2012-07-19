function yi = quick_interp(x, y, xi,varargin)
% yi = quick_interp(x, y, xi)
%
% Quick interpolation. Assumes x is monotonically increasing. Returns
% first/last y-value for out-of-range xi.

% if xi>x(end)
%     yi = y(end);
%     return
% elseif xi<x(1)
%     yi = y(1);
%     return
% else
%     yi = lininterp1f(x,y,xi,0);
%     return
% end

%  Legacy method, use this if lininterp1f fails.
ind = find(x>=xi,1);

% if isempty(ind)
%     yi = y(end);
%     return
% elseif ind == 1
%     yi = y(1);
%     return
% else
%     yi = (y(ind)-y(ind-1))*(xi-x(ind-1))/(x(ind)-x(ind-1))+y(ind-1);
% end

if isempty(ind) || ind == 1
    if isempty(varargin)
        if isempty(ind)
            yi = y(end);
            return
        elseif ind == 1
            yi = y(1);
            return
        end
    else
        p = polyfit(x,y,1);
        yi = polyval(p,xi);
        return
    end
else
    yi = (y(ind)-y(ind-1))*(xi-x(ind-1))/(x(ind)-x(ind-1))+y(ind-1);
end