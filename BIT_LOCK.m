function index = BIT_LOCK(phi, prn)
% index = BIT_LOCK(phi, prn)
% 
% Finds the index of the first chip in the first whole bit, using the
% histogram method.

bit = (phi > pi/2) | (phi < -pi/2);

transition = xor(bit(1:end-1),bit(2:end));

L = length(transition);
w = floor(L/20);

matrix = reshape(transition(1:w*20),20,w);

histogram = sum(matrix,2);

[N, I] = max(histogram);

% It is nice to get a warning if a bit lock is unreliable. This model
% assumes that all peaks are binomial, except the highest. The parameter p
% is estimated from all but the highest peak. The probability that the
% maximum of 20 peaks would be higher than the actual height of the
% heighest peak is then calculated.

% Estimate P.
P = sum(histogram((1:20)~=I))/(19*L);

% Calculate P[Max(Bin(L;P)_20)>N].
pNoise = 1-binocdf(N,L,P)^20;

% If there is a high risk of noise, warn the user.
if pNoise>0.01
    warning('There is a %.0f%% risk that the BST is corrupted by noise for SV%02i.',...
        pNoise*100,prn)
end

% Convert to Matlab index.
index = mod(I,20)+1;