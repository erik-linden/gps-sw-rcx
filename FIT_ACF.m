function [cst_rel, phase, magnitude, lock, rmse] = FIT_ACF(IQ, type)
% [cst_rel, phase, magnitude, lock, rmse] = FIT_ACF(IQ, type)
%
% Measures cst using curve fitting.
%
% Input:
%     IQm        magnitude of the output from CORRELATE
%
% Output:
%     cst_rel     cst relative the center correlator
%     lock        some sort of lock information

global N_CORR MAX_SPACE CHIP_RATE

index = -(N_CORR-1)/2:(N_CORR-1)/2;     %index for correlators

switch type
    case 'tent'
        % TENT FUNCTION
        
        [IQmax indMax] = max(abs(IQ));  %used in setting the start point
        phaseZero = angle(IQ(indMax));
        
        [magnitude, phase, center, lock, rmse] = fit_tent(index',IQ,IQmax,phaseZero,index(indMax),N_CORR/MAX_SPACE);
        
        if lock
            cst_rel = center*MAX_SPACE/(CHIP_RATE*(N_CORR-1));  %convert from index to time
        else
            cst_rel = NaN;
        end
        
    case 'poly'
        % QUADRATIC FUNCTION
        
        p = polyfit(index',abs(IQ),2);
        
        if p(1) >= 0
            % The arc is curving upwards.
            
            lock = false;
            cst_rel = NaN;
            magnitude = NaN;
            rmse = NaN;
            
        elseif -p(2)/(2*p(1)) < index(1) || -p(2)/(2*p(1)) > index(end)
            % The center is not in the range.
            
            lock = false;
            cst_rel = NaN;
            magnitude = NaN;
            rmse = NaN;
            
        else
            % Otherwise, the lock is valid.
            
            lock = true;
            center = -p(2)/(2*p(1));
            cst_rel = center*MAX_SPACE/(CHIP_RATE*(N_CORR-1));
            maximum = (-3/4*p(2)^2/p(1)+p(3));
            rmse = sqrt(mean((IQm-polyval(p,index')).^2));
            
            % The conversion factor here is the get the actual peak hight.
            % It is based on analytical calculations of a quadratic
            % function fitted to a tent function. See the relevant
            % Mathematica notebook.
            magnitude = maximum*32/(32-3*MAX_SPACE);
            
            % This is all to interpolate the phase.
            n_max = cst_rel*(CHIP_RATE*(N_CORR-1)/MAX_SPACE)+(N_CORR+1)/2; %get the cst in "index time"
            n_f = floor(n_max);   	%get samples to the left and right of the cst
            n_c = ceil(n_max);
            
            phi = unwrap(angle(IQ([n_f n_c])));  %get phase at samples
            if n_f ~= n_c
                phase = (phi(2)-phi(1))*(n_max-n_f)/(n_c-n_f)+phi(1);
            else
                phase = phi(1);
            end
        end
        
    otherwise
        error('No curve fitting function "%s"', type)
end
