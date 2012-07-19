function uvw = project_vel(xyz, pos)
% uvw = project_vel(xyz, pos)
% 
% Projects vectors xyz onto the local Earht plane, at pos. Differs from
% project_to_surface in that the w-vector is not treated as altitude.

dim = size(xyz,1);

% Set up basis vectors.
% X+ -> East
% Y+ -> North
% Z+ -> Up
b1 = [-pos(2)/pos(1),1,0];
b1 = b1./norm(b1);
b2 = -[pos(3),pos(2)*pos(3)/pos(1),-pos(1)-pos(2)^2/pos(1)];
b2 = b2./norm(b2);
b3 = [pos(1),pos(2),pos(3)];
b3 = b3./norm(b3);

uvw = zeros(dim,3);

for n = 1:dim
    uvw(n,1) = dot(b1,xyz(n,1:3));
    uvw(n,2) = dot(b2,xyz(n,1:3));
    uvw(n,3) = dot(b3,xyz(n,1:3));
end