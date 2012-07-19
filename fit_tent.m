function [magnitude, phase, center, lock, rmse] = fit_tent(x,y_c,magnitude,phase,center,width)
% [magnitude, center, lock, rmse] = fit_tent(x,y,magnitude,center,width)

global COH_INT_TIME SEARCH_WIN

% persistent shape
% if isempty(shape)
%     shape = zeros(2*64+1,1);
% end

% These parameters control the fitting.
maxEval = 12;                   % Typically 3-5, but much larger when there are convergence problems
rangeErr = 1e-3;                % Assume convergence if change in cst is less than this, in meters

% Initialize values.
a = magnitude;                  % Get the initial hight of the "tent"
c = center;                     % Get the initial center
k = width;                      % Get the initial width
p = phase;                      % Get the initial width
W = diag(repmat(abs(x-c)<SEARCH_WIN*width/2,[2,1]));

% These control when to abort the fitting.
aMin = 0.67*noise_std(COH_INT_TIME);     % Break if the signal falls below the noise 50% limit
minDeltaC = width*rangeErr/(3e8/1023e3); % This is the minimun cst change, in samples

% Initialize loop variables.
nPoints = length(x);            % Number of data points
y = zeros(2*nPoints,1);
y(1:nPoints) = real(y_c);
y(nPoints+1:end) = imag(y_c);
beta = [magnitude, center, phase];

% 'TolX', minDeltaC,
options = optimset(...
    'Display','off',...
    'MaxFunEvals',25,...
    'Jacobian','on',...
    'Algorithm','levenberg-marquardt');
[beta,resnorm,~,flag] = lsqnonlin(@getRes,beta,...
    [],[],options);

a = beta(1);
c = beta(2);
p = beta(3);
% 
% % ## PLOTTING ##
% % Used to monitor fitting.
% plot(x/width,real(exp(-1i*p)*y_c)/a)
% hold on
% plot(x/width,imag(exp(-1i*p)*y_c)/a,'r')
% plot(x/width,tent(x, c, 1, k),'g')
% % plot(c,1,'.b')
% ylim([-.1, 1])
% hold off
% ylim([-0.1,1.1])
% xlim([-1.5,1.5])
% legend({'I-channel','Q-channel','Curve fit'})
% xlabel('Delay - chips')
% ylabel('Relative amplitude')
% drawnow
% pause(1/24)
% 
% % Used to monitor signal "shape", uncomment persistent on first
% % lines.
% shape = [shape,interp1(x-center,exp(-1i*p)*y_c,x,'nearest',0)];
% perI = prctile(real(shape),[50-34.1 50 50+34.1],2);
% scl = max(perI(:,2));
% plot(x/width,perI(:,2)/scl,'b')
% hold on
% perQ = prctile(imag(shape),[50-34.1 50 50+34.1],2);
% plot(x/width,perQ(:,2)/scl,'r')
% 
% 
% plot(x/width,perI(:,1)/scl,'b:')
% plot(x/width,perI(:,3)/scl,'b:')
% 
% plot(x/width,perQ(:,1)/scl,'r:')
% plot(x/width,perQ(:,3)/scl,'r:')

% shape = shape + interp1(x-center,exp(-1i*p)*y_c,x,'nearest',0);
% plot(x/width,real(shape)/max(real(shape)))
% hold on
% plot(x/width,imag(shape)/max(real(shape)),'r')
% hold off
% ylim([-0.1,1.1])
% xlim([-1.5,1.5])
% legend({'I-channel','Q-channel'})
% xlabel('Delay - chips')
% ylabel('Relative amplitude')
% drawnow;

% Decide what to return.
aGood = a>aMin;
% cGood = min(x-c)<-SEARCH_WIN*width/2 && ...
%     max(x-c)>SEARCH_WIN*width/2;
cGood = min(c)>min(x) && max(c)<max(x);
flagGood = flag>=0;

if aGood && cGood && flagGood
    magnitude = a;
    phase = p;
    center = c;
    lock = true;
    rmse = sqrt(sqrt(resnorm/nPoints)); %RMS of the _complex_ error vector
else
    % If we fell through the loop or it was broken, we have no lock.
    phase = NaN;
    lock = false;
    rmse = NaN;
    magnitude = NaN;
    center = NaN;
end

return

    function [R,J] = getRes(beta)
        a = beta(1);
        c = beta(2);
        p = beta(3);
        
        y_mr = cos(p)*tent(x, c, a, k);
        y_mi = sin(p)*tent(x, c, a, k);
        y_m = [y_mr; y_mi];
        R = y_m - y;
        
        x_tent = (x>c-k) & (x<c+k);
        x_neg = x<c;
        
        J = zeros(2*nPoints,3);
        J(1:nPoints,1)     =  cos(p).*x_tent.*(1-1/k*(x-c).*(-1).^x_neg);     %dyda-I
        J(1:nPoints,2)     =  cos(p).*x_tent.*a./k.*(-1).^x_neg;              %dydc-I
        J(1:nPoints,3)     = -sin(p).*x_tent.*a.*(1-1/k*(x-c).*(-1).^x_neg);  %dydp-I
        J(nPoints+1:end,1) =  sin(p).*x_tent.*(1-1/k*(x-c).*(-1).^x_neg);     %dyda-Q
        J(nPoints+1:end,2) =  sin(p).*x_tent.*a./k.*(-1).^x_neg;              %dydc-Q
        J(nPoints+1:end,3) =  cos(p).*x_tent.*a.*(1-1/k*(x-c).*(-1).^x_neg);  %dydp-Q
        
        R = W*R;
        J = W*J;
    end
end