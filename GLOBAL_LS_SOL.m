
sv = [2 5 10 12 21 23 25 29 30 31];%
% sv = [5 7 8 9 11 15 17 18 19 26 27 28];%
% sv = [2 5 6 10 12 13 16 21 23 25 29 30 31];
% sv = [3 5 7 8 15 18 19 21 26 28];
% sv = [2 5 10 12 21 23 25 29 30 31];
r = 1;

w = [1/10,1/0.2,1/1,1e-2].^2;

[x time]=nav_residuals('Init',sv,100,[w(1:3),w(4)],[0,0,0]);
lb = -inf(size(x));
lb(end-5:end-3) = 0.9;
lb(end-2:end)   = -1;

ub = +inf(size(x));
ub(end-5:end-3) = 1.1;
ub(end-2:end)   = +1;

options = optimset('Jacobian','on',...
    'TolX',1e-8,...
    'TolFun',1e-8,...
    'MaxIter',160);
x = lsqnonlin(@(x) nav_residuals('Func',x),x,lb,ub,options);
% nav_residuals('turnOn',[1,0,0]);
% [x,~,residual] = lsqnonlin(@(x) nav_residuals('Func',x),x,lb,ub,options);
nav_residuals('turnOn',[1,0,0]);
[x,~,residual] = lsqnonlin(@(x) nav_residuals('Func',x),x,lb,ub,options);
% nav_residuals('turnOn',[1,1,1]);
% x = lsqnonlin(@(x) nav_residuals('Func',x),x,lb,ub,options);
nav_residuals('Func',x)


[x2 time2]=nav_residuals('Init',sv,r,w,[1,0,0]);

N1 = length(time);
N2 = length(time2);

for i = 1:4
    x2(N2*(i-1)+(1:N2)) = interp1(time,x(N1*(i-1)+(1:N1)),time2,'linear','extrap');
end
 x2(N2*4+6+(1:N2-1)) = interp1(time(1:end-1),x(N1*4+6+(1:N1-1)),time2(1:end-1),'linear','extrap');

lb = -inf(size(x2));
lb(4*N2+(1:3)) = 0.9;
lb(4*N2+(4:6)) = -3;

ub = +inf(size(x2));
ub(4*N2+(1:3)) = 1.1;
ub(4*N2+(4:6)) = +3;

[x2,~,residual] = lsqnonlin(@(x) nav_residuals('Func',x),x2,lb,ub,options);
% nav_residuals('turnOn',[1,1,1]);
% [x2,~,residual] = lsqnonlin(@(x) nav_residuals('Func',x),x2,lb,ub,options);
nav_residuals('Func',x2)