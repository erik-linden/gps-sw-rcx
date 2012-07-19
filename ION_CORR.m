function [tDelay, tDcb] = ION_CORR(gpsTime,recPos,svPos,prn)

% This is used to find the inteception point.
earthRadius = 6371009;

global INITIALIZE RECEIVER_FILE
persistent ION

% Load TEC data.
if ~exist('ION','var') || INITIALIZE
    load(RECEIVER_FILE);
    if Rec.IonType == 0
        error('RAIN:NavData','Ion LUT is empty.')
    end
    ION = Rec.Ion;
    clear('Rec');
end

% Find the first map after the current time.
ind = find(ION.Time>gpsTime,1);

% Read the maps before and after the current time.
map1 = ION.TEC(:,:,ind-1);
map2 = ION.TEC(:,:,ind);

% Compute the sideral time for the current time and for the two maps.
t = gpsTime/240;
Ti = ION.Time(ind-1)/240;
Tip = ION.Time(ind)/240;

% Radius of the interception point.
rIP = ION.Height + earthRadius;

% Radius of the receiver.
rREC = sqrt(sum(recPos.^2));

% Compute the location of the intercept-point.
d = svPos-recPos;
dsqr = (norm(d)^2);
x = dot(d,recPos)/dsqr;

vecIP = recPos + (sqrt(x^2+(rIP^2-rREC^2)/dsqr)-x)*d;

% This is the cosine of the delflection of the receiver/IP vector.
cosDef = dot(recPos,d)/(norm(recPos)*norm(d));

% Weight the TEC based on the grazing angle.
weight = (1 - (sqrt(1-cosDef^2)*rREC/rIP)^2)^(-1/2);

% Get the IP in lon/lat.
latlonalt = ecef2lla(vecIP);
lat = latlonalt(1);
lon = latlonalt(2);

% Get the longitude of the IP in the two maps, using a sun-fixed coordinate
% system.
lp = lon+(t-Ti);
lpp = lon+(t-Tip);

% Interpolate the TEC value at the IP location, in both maps.
[X,Y] = meshgrid(ION.LatGrid,ION.LonGrid);

Ei = interp2(X,Y,map1,lat,lp);
Eip = interp2(X,Y,map2,lat,lpp);

% Interpolate the TEC value in time.
E = (Tip-t)/(Tip-Ti)*Ei+(t-Ti)/(Tip-Ti)*Eip;

% This is the L1 group delay.
rhoL1 = weight*1.34e9*E/((10.23e6*154)^2);

% Read the DCB from the ion data.
DCB = -1.55*ION.DCB(prn);

% Output the results.
tDelay = rhoL1;
tDcb = DCB;