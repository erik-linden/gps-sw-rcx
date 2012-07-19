function y = tent(x, center, height, half_width)
% y = tent(x, center, height, half_width)
%
% Returns a tent function, __/\__

x_tent = (x>center-half_width) & (x<center+half_width);
x_neg = x<center;
y = x_tent.*height.*(1-1/half_width.*(x-center).*(-1).^x_neg);

