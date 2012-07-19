function [IQ t_center] = CORRELATE(in_sig, cst, w_if, phi_if, corr, indexRelSOB)
% [IQ t_center] = CORRELATE(in_sig, cst, w_if, phi_if, corr)
%
% Correlates a given signal with a set of code replicas.
%
% Input:
%     in_sig      - signal to search for code in
%     cst         - estimated code start time, relative beginning of in_sig
%     w_if        - estimated angular intermediate frequency at cst
%     phi_if      - estimated intermediate frequency phase at cst
%     corr        - matrix of code replicas [COH_INT_SAMxN_CORR]
%
% Output:
%     IQ          - [N_CORRx1] complex correlator output for each code replica
%     t_center    - the time which the correlations are done relative to
%
% We can only do correlations at actual samples, and the cst will generally
% not coincide with a sample. So we find the closest sample and then output
% the cst we actually used as t_center, relative the start of the file.

global INITIALIZE COH_INT_SAM COH_INT_TIME ONE_MSEC_SAM TP

persistent polarityUnknown isInv deltaMax

% Initialize persistent variables.
if isempty(polarityUnknown) || INITIALIZE
    % Assume unknown polarity to start with.
    polarityUnknown = true;
    
    % Compute how much larger one polarity must be than the other to be
    % able to separat them with a given degree of statistical confidence.
    deltaMax = sqrt(1-pi/4)*noise_std(COH_INT_TIME)*norminv(.9995,0,1);
end

I = round(cst/TP);  %find sample nearest to cst

t_center = I*TP;    %report the actual cst used

samp_index = I:I+COH_INT_SAM-1;    %prepair a index vector

time = (samp_index*TP-cst)';       %make a time vector. time(0) is not necissarly 0, so phase will be correct.
freq_arg = time*w_if + phi_if;     %get frequency argument

ii = 1:COH_INT_SAM; %to use for the correlation

I = I+1;            %matlab indexes from 1

BB = in_sig(I+ii-1).*(cos(freq_arg)-1j*sin(freq_arg));  %mix to complex baseband

% This part does the correlations and uses a lot of logic to handle bit
% flips. It needs to decide if the correlator might be straddling a bit
% flip an, in that case, see if the polarity of the flip already has been
% resolved or if it needs to be tested.

% This is the number of samples left before the possible bitflip.
samplesBeforeInversion = (20-indexRelSOB)*ONE_MSEC_SAM;

% If the number of samples in the correlator is more than the number of
% samples left, we are straddling a bit flip.
isStraddlingBitFlip = samplesBeforeInversion<COH_INT_SAM;

% If we are not straddling or know that the bit don't flip or the flip is
% unknown, use the  un-fliped correlator.
if ~isStraddlingBitFlip || (isStraddlingBitFlip && (polarityUnknown || ~isInv))
    IQ = corr*BB;
end

% If we are straddling and know the bit is inverted or if the polarity of
% the flip is unknown, use the flipped correlator.
if isStraddlingBitFlip && (polarityUnknown || isInv)
    prt1 = corr(:,1:samplesBeforeInversion)*BB(1:samplesBeforeInversion);
    prt2 = corr(:,samplesBeforeInversion+1:end)*BB(samplesBeforeInversion+1:end);
    IQinv = prt1-prt2;
end

% If we are not straddling, polarity of the next flip is unknown.
if ~isStraddlingBitFlip
    polarityUnknown = true;
    return
    % If we are straddling and the bit isen't inverted, return uninverted
    % correlator.
elseif ~polarityUnknown && ~isInv
    return
    % If we are straddling and the bit is inverted, return inverted
    % correlator.
elseif ~polarityUnknown && isInv
    IQ = IQinv;
    return
    % If the polarity is unknown, test by looking at the maximum.
else
    peakMagStr = max(abs(IQ(:)));
    peakMagInv = max(abs(IQinv(:)));
    
    % Decide if polarity can be determined.
    if abs(peakMagInv-peakMagStr)>deltaMax
        polarityUnknown = false;
    end
    
    % Take the largest polarity as the result.
    if peakMagInv>peakMagStr
        isInv = true;
        IQ = IQinv;
    else
        isInv = false;
    end
    
    return
end
