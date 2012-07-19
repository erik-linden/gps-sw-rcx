function [velocity, lo, dop, res] = SOLVE_VEL(sv, cpr, recPos, gpsTime, weights)
% [velocity, lo, dop, res] = SOLVE_VEL(sv, cpr, recPos, gpsTime, weights)

global C L1

% Initialize SV-related vectors.
nSv = length(sv);
svPos = zeros(nSv,3);
svVel = zeros(nSv,3);
svCarCorr = zeros(nSv,1);

for n = 1:nSv
    % Get aprox. position.
    svPosTemp = GET_EPH(gpsTime,sv(n),'pos');
    
    % Get the range.
    range = sqrt(sum((svPosTemp-recPos).^2,2));
    
    % Compute the time of transmission;
    gpsTimeAtTx = gpsTime-range/C;
    
    % Do a total of three fixed-point iterations.
    for k = 1:2
        svPosTemp = GET_EPH(gpsTimeAtTx,sv(n),'pos');
        range = sqrt(sum((svPosTemp-recPos).^2,2));
        gpsTimeAtTx = gpsTime-range/C;
    end       
    
    % Get position and velocity.
    svPos(n,:)   = GET_EPH(gpsTimeAtTx,sv(n),'pos');
    svVel(n,:)   = GET_EPH(gpsTimeAtTx,sv(n),'vel');
    
    % We also want to correct for any error in the L1 carrier. I simply
    % assume that the error is the time derivative of the clock error,
    % times the nominal carrier frequency.
    svCarCorr(n) = GET_EPH(gpsTimeAtTx,sv(n),'dClkdt')*L1;
end

% Initialize velocity vectors.
uVec = nan(nSv,3);
vel_sv_rad = nan(nSv,1);
vel_rec_rel = nan(nSv,1);

for n = 1:nSv
    % Relatevistic gamma
    gamma = 1/sqrt(1-(norm(svVel)/C)^2);
    
    % Create a unit vector from receiver to SV
    uVec(n,:) = recPos - svPos(n,:);
    uVec(n,:) = uVec(n,:)./norm(uVec(n,:));
    
    % Project the velocities onto the generated basis vector (positiv =>
    % separating).
    vel_sv_rad(n) =  dot(uVec(n,:),svVel(n,:))/gamma;
    
    % Finally, compute doppler frequency.
    vel_rec_rel(n) = -(cpr(n)/(2*pi) - svCarCorr(n))*C/L1;
end

% Create relative velocity vector
l = vel_sv_rad+vel_rec_rel;

% Create matrix E.
E = [uVec, ones(nSv,1)];

% Solve the system.
V = lscov(E,l,weights);

% Compute velocity residuals.
res = E*V-l;

% Compute (weighted) DOP.
dop = diag(inv(E'*diag(weights)*E));

% Exctract velocity and LO-offset.
velocity=V(1:3);
lo=V(4)*(2*pi*L1)/C;