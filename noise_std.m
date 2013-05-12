function s = noise_std(T)
%s = noise_std(T)
%
% Calculate noise sigma based on an sqrt power law, sigma = a*T^(1/2)

global INITIALIZE RECEIVER_FILE
persistent alpha

%Make sure LUT is loaded.
if ~exist('alpha','var') || INITIALIZE
    if(~isempty(RECEIVER_FILE))
        load(RECEIVER_FILE);
        alpha = Rec.Alpha;
        clear('Rec');
    else
        alpha = 5e3;
    end
end

s = alpha*T^0.5;
