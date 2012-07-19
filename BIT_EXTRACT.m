function bit = BIT_EXTRACT(phi, I)
% bit = BIT_EXTRACT(phi, I)
% 
% Converts chips to bits by taking the mean.

dim1 = floor(length(phi)/20)-2;

bit = zeros(dim1, 1);

for m = 1:dim1
    index = I+((m-1)*20:m*20-1);
    
    bit(m) = nanmean(mod(phi(index)+pi/2,2*pi))>pi;
end
