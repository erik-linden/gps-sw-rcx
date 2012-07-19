function Coef = GEN_EPH_LUT(files)
% Coef = GEN_EPH_LUT(files)
%
% Generates LUT for ephemerideis data. Parameters should not be changed.
% The coefficients fits a N:degree polynom to N+1 data points
%
% Input
%         files:  {'file1', 'file2'...}
%
% Output
% Coef =
%       GPSTime:  [1xEpocs]         %system time at the middle of segments
%             X:  [32xEpocsxCoef]   %coeffiecients
%             Y:  [32xEpocsxCoef]
%             Z:  [32xEpocsxCoef]
%          dClk:  [32xEpocsxCoef]
%             N:  [1x1]             %interpolation parameters
%             M:  [1x1]

% Parameters
N = 10;  %Polynom degree (16), must be even
M = 4;    %Points in valid region (7), must be odd
Coef = {};

S = READ_SP3(files);%Get a struct of ephemeridies data

% Quit if S is empty.
if isempty(S)
    warning('RAIN:NavData','SP3 struct is empty.');
    return
end

x = (-N/2:1:N/2);   %Generate sample points
x = 2*x/N;          %Scale x to [-1 1], this improves numerics

V = vander(x);      %Generate vandermond matrix

time = S.GPSTime;   %vector of times

segments = floor((length(time)-(N+1))/(M-1)); %Number of segments

for n = 0:segments
    for prn = 1:32
        Coef.GPSTime(n+1) = time((M-1)*n+N/2+1); %time at the middle of the segment
        index = (M-1)*n+(1:N+1);    %index of points to interpolate
        
        Coef.X(prn,n+1,1:N+1)    = V\S.X(prn,index)';
        Coef.Y(prn,n+1,1:N+1)    = V\S.Y(prn,index)';
        Coef.Z(prn,n+1,1:N+1)    = V\S.Z(prn,index)';
        Coef.dClk(prn,n+1,1:N+1) = V\S.dClk(prn,index)';
    end
end

Coef.N = N; %Store parameters in the LUT
Coef.M = M;