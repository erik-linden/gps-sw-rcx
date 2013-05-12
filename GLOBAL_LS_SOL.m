function GLOBAL_LS_SOL(sv,settings)

% Standard deviations for different residuals.
rangeStd    = 10;  % Range std in meters
velStd      = 0.2; % Doppler std in m/s
accStd      = 1.0; % Acceleration magnitude std in m/s^2
jerkStd     = 7e2; % Jerk std in m/s^3
cJerkStd    = inf; % Clock jerk in m/s^3
wuaStd      = 0.2/5;% Wind-up acc. in m/s^2

% Set lsqnonlin settings.
options = optimset('Jacobian','on',...
        'TolX',1e-8,...
        'TolFun',1e-8,...
        'MaxIter',160,...
        'Display','notify');

% Give the user a hint about the weighting scheme.
fprintf(['\nStandard deviations are:\nRange \t\t%4.2g m\nVel \t\t%4.2g m/s\n'...
    'Acc \t\t%4.2g m/s^2\nJerk \t\t%4.2g m/s^3\nClock jerk \t%4.2g m/s^3\nWind-up acc %4.2g m/s^2\n'],...
    rangeStd,velStd,accStd,jerkStd,cJerkStd,wuaStd)
disp('These values are hard coded in GLOBAL_LS_SOL')
disp(' and there is no reason to assume that they are good.')

nSteps = size(settings,1);
fprintf('\nSolving in %u steps.\n',nSteps)

% Weight vector, inverse varaince.
w = 1./[rangeStd, velStd, accStd, jerkStd, cJerkStd, wuaStd].^2;
x_old = [];

for n = 1:nSteps
    r = settings(n,1);
    enable = settings(n,2:4);
    
    fprintf('\nInitializing at %u Hz with switch [Vel:%u Acc:%u Jerk:%u]\n',1000/r,enable)
    [x time]=nav_residuals('Init',sv,r,w,enable);
    N = length(time);
    
    % Upper and lower bounds for acc. bias and slope.
    lb = -inf(size(x));
    lb(4*N+(1:3)) = 0.95;
    lb(4*N+(4:6)) = -1;

    ub = +inf(size(x));
    ub(4*N+(1:3)) = 1.05;
    ub(4*N+(4:6)) = +1;
    
    % Interpolate from old solution, if avialable.
    if(~isempty(x_old))
        fprintf('\nInterpolating old solution.\n')
        
        N_old = length(time_old);
        
        for i = 1:4
            % Interpolate the 3 position varaibles and the one clock error variable
            x(N*(i-1)+(1:N)) = ...
                interp1(time_old,x_old(N_old*(i-1)+(1:N_old)),time,'linear','extrap');
        end
        % Interpolate 6 acc. bias variables.
        x(N*4+(1:6)) = x_old(N_old*4+(1:6));
        % Interpolate the wind-up.
        x(N*4+6+(1:N-1)) = ...
            interp1(time_old(1:end-1),x_old(N_old*4+6+(1:N_old-1)),time(1:end-1),'linear','extrap');
        
    end
    
    fprintf('\nSolver started.\n')
    x = lsqnonlin(@(x) nav_residuals('Func',x),x,lb,ub,options);
    
    % Store new solution.
    x_old = x;
    time_old = time;
    
end

% Save solution.
fprintf('\nSaving solution\n')
nav_residuals('Func',x)

% Save used settings.
save(sprintf('%sposLS_settings',TRACK_DIRECTORY),'rangeStd',...
    'velStd','accStd','jerkStd','cJerkStd','wuaStd')

