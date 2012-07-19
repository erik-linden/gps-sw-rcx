function [lock, rmse, magnitude, cst_actual, cst_err, phi_actual, phi_err] = MEASURE(...
    in_sig, cst, w_if, phi_if, corr, t0, indexRelSOB)
% [lock, magnitude, cst_actual, cst_err, phi_actual, phi_err] = MEASURE(in_sig,
% cst, w_if, phi_if, corr, t0)
%
% Takes predictions of cst, w_if and phi_if and returns measurments of cst
% and phi_if, with the associated errors.
%
% Input:
%     in_sig      - vector of recorded data
%     cst         - predicted cst
%     w_if        - predicted angular intermediate frequency at predicted cst
%     phi_if      - predicted intermediate frequency phase at predicted cst
%     corr        - [COH_INT_SAMxN_CORR] matrix of upsampled correlators, from GEN_CORR
%     t0          - start time of current file, for seemless multiple file access
%
% Output:
%     cst_actual  - measured cst
%     cst_err     - difference between predicted and measured cst
%     phi_actual  - measured phase at measured cst
%     phi_err     - difference between predicted and measured phase
%     magnitude   - an estimate of signal magnitude
%     lock        - curve fitting worked/didn't work

global CF_TYPE

cst_file = cst - t0; %get time relative start of file

% Call the function that does the correlations.
[IQ t_center] = CORRELATE(in_sig, cst_file, w_if, phi_if, corr, indexRelSOB);

% Call the function that does the curve fitting on the output from the
% correlator.
[cst_rel, phi_actual, magnitude, lock, rmse] = FIT_ACF(IQ, CF_TYPE);

if lock
    % If curve fitting locked, return data...

    cst_actual = t0 + t_center + cst_rel;   % calculate the actual cst...
    cst_err = cst - cst_actual;             % ...and the error

    phi_err = mod(phi_actual+pi/2,pi)-pi/2; % this is equivalent to an arctan discriminator
else
    % If not, return NaNs.

    cst_actual  = NaN;
    cst_err     = NaN;
    phi_actual  = NaN;
    phi_err     = NaN;
    magnitude   = NaN;
end