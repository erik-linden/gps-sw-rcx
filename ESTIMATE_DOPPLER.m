function [df, varargout] = ESTIMATE_DOPPLER(t, prn)
% w_df = ESTIMATE_DOPPLER(t, prn)
%
% Returns expected doppler frequency based on ephmeries LUT
% and receiver data file.
% 
% Input
% 	t       - time relative start of recorded data
% 	prn     - SV to get doppler for
%
% Output
%   df      - derived doppler frequency

L1 = (10.23e6*154); %L1 carrier, Hz
C = 299792458; %speed of light

% Get GPS system time from receiver data.
gps_time = GET_REC_DATA(t, 'time');

% Get receiver postion and velocity.
recXYZ = GET_REC_DATA(t, 'pos');
recUVW = GET_REC_DATA(t, 'vel');

% Get SV postion and velocity.
svXYZ = GET_EPH(gps_time, prn, 'pos');
svUVW = GET_EPH(gps_time, prn, 'vel');

% Create a unit vector from receiver to SV
uVec = svXYZ - recXYZ;
uVec = uVec./norm(uVec);

% Project the velocities onto the generated basis vector (positiv =>
% separating).
vel_rel =  dot(uVec,svUVW) - dot(uVec,recUVW);

% Finally, compute doppler frequency.
df = -vel_rel*L1/C;

if nargout == 2
    elev = 90-180/pi*(acos(dot(uVec,recXYZ)/(norm(recXYZ))));
    varargout = {elev};
end

