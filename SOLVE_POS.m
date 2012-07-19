function [posObs, gpsTimeOffset, dop, res] = SOLVE_POS(sv,rawRange,svTime,weights,guess,useIon)
% [posObs, gpsTime, dop, res] = SOLVE_POS(sv,rawRange,svTime,weights,guess)

global C

% Speed of light and earth rotation speed.
OmegaE = 7.2921151467e-5;

% Initialize vectors.
nSv = length(sv);

corrRange = zeros(nSv,1);
range = zeros(nSv,1);
satXYZ = zeros(nSv,3);
gpsTime = zeros(nSv,1);
guessmatrix = zeros(nSv, 3);
cmatrix = ones(nSv,1);

% Unset stop flag.
stop = false;

% Initilalize correction vectors.
tRel=zeros(nSv,1);
tIon=zeros(nSv,1);
tDcb=zeros(nSv,1);

% Load APC corrections.
apc = APC_CORR(sv);

% Get corrections.
for n = 1:nSv
    
    % The clock error correction is given in gps-time, as
    %  T_SV - dClk(T_GPS) = T_GPS
    %
    % This is solved by fixed point iteration.
    dClk = 0;
    for k = 1:2
        dClk = GET_EPH(svTime(n) - dClk,sv(n),'dClk');
    end
    
    % Compute the GPS-time for the current SV, at the time of transmition.
    gpsTime(n) = svTime(n) - dClk;
    
    % Compute the position and velocity of the current SV at the time of
    % transmition.
    pos = GET_EPH(gpsTime(n),sv(n),'pos');
    vel = GET_EPH(gpsTime(n),sv(n),'vel');
    
    % Compute the relativistic correction.
    tRel(n) = -2*dot(pos,vel)/C^2;
    
    if useIon
        % Get the ion delay and differential-code-bias.
        [t1 t2] = ION_CORR(svTime(1),guess,pos,sv(n));
        tIon(n) = t1;
        tDcb(n) = t2;
    else
        tIon(n) = 0;
        tDcb(n) = 0;
    end
    
    % Compute the total time correction.
    tCorr = dClk + tRel(n) + tIon(n) - tDcb(n);
    
    % Compute the range correction, including APC correction.
    corrRange(n) = (rawRange(n)' + tCorr)*C + apc(n);
end

while (stop == 0)
    % create guess matrix
    for i = 1:nSv
        guessmatrix(i,:) = guess;
    end
	% calculate difference between time of reception and
    % transmission in order to shift satellites back to where
    % they were at transmission
    tcorr = range./C;
    % calculate 'satXYZ' based at time 'GPStime' - 'tcorr' for
    % each satellite; then, rotate backwards by Earth's rotation
    % during that time
    for i = 1:nSv
        satXYZ(i,:) = GET_EPH(gpsTime(i),sv(i),'pos');
    end
    satX = satXYZ(:,1);
    satY = satXYZ(:,2);
	satZ = satXYZ(:,3);
    delX =  cos(tcorr*OmegaE).*satX + sin(tcorr*OmegaE).*satY;
    delY = -sin(tcorr*OmegaE).*satX + cos(tcorr*OmegaE).*satY;
    satXYZ = [delX delY satZ];

    % calculate satellite ranges 'range' from guess to satellites
    delXYZ = satXYZ - guessmatrix;
    range = sqrt(sum((delXYZ.^2),2));
    % form the vector 'l' and the matrix 'A'
    l = corrRange - range;
    pmatrix = range*ones(1,3);
    A = delXYZ ./ pmatrix;
    A = -[ A cmatrix ];
    % get weights
    % solve for  'deltaPos' which is contains dx, dy, dz, and  dt
    deltaPos = lscov(A,l,weights);
    % calculate 'obsPos' by adding 'deltaPos' to the current guess
    posObs = [guess 0] + deltaPos';
    posObs(4) = posObs(4)/C;
    % check to see if the initial guess and the computed result is
    % "close enough"; if it is, then stop the iteration by setting the
    % 'stop' flag to 1; else, set the new initial guess equal to the
    % last computed value, and iterate again
    if (norm(posObs(1:3) - guess) < 1e-6)
        stop = 1;
    end
    guess = posObs(1:3);
end
% GPS-time...
gpsTimeOffset = posObs(4);
% Position...
posObs = posObs(1:3);
% Equation residuals...
res =A*deltaPos-l;
% Weighted dilution-of-precision...
dop = diag(inv(A'*diag(weights)*A));