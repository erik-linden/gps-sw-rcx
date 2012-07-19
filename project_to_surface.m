function [uvw] = project_to_surface(xyz, varargin)
% [uvw] = project_to_surface(xyz, [center])
% 
% Projects the ECEF vector xyz onto the local Earth surface plane.

if ~isempty(varargin)
    pos = varargin{1};
else
    pos = nanmean(xyz);
end

dim = size(xyz,1);

% Set up basis vectors.
% X+ -> East
% Y+ -> North
% Z+ -> Up
b1 = [-pos(2)/pos(1),1,0];
b1 = b1./norm(b1);
b2 = -[pos(3),pos(2)*pos(3)/pos(1),-pos(1)-pos(2)^2/pos(1)];
b2 = b2./norm(b2);
% b3 = [pos(1),pos(2),pos(3)];
% b3 = b3./norm(b3);

uvw = zeros(dim,3);

for n = 1:dim
    uvw(n,1) = dot(b1,xyz(n,1:3));
    uvw(n,2) = dot(b2,xyz(n,1:3));
    %uvw(n,3) = dot(b3,xyz(n,1:3)); % Not used, altitude instead
end

lla = ecef2lla(xyz(:,1:3));
uvw(:,3) = lla(1:dim,3);