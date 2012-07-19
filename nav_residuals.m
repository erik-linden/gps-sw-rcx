function  varargout = nav_residuals(string,varargin)

persistent rangeCorr wRange recTime prRefTime accBody accCorr
persistent dopplerCorr wDoppler svPos svVel rNacc
persistent deltaT N nSv jRes jPar jVal rNtot rNwuv
persistent rNrange rNvel rNjerk  w sqrW svPosCorr svVelCorr

switch string
    case 'Init'
        % Read out the content of the input arguments.
        sv = varargin{1};
        r  = varargin{2};
        weights = varargin{3};
        turnOn  = varargin{4};
        
        % Load all the data.
        [rangeCorr, wRange, recTime, prRefTime,...
            dopplerCorr, wDoppler, svPos, svVel] = CORRECTED_MES(sv,r);
        
        [accBody accCorr] = GET_ACCELERATION(recTime);
        
        nSv = length(sv);
        N = length(recTime);
        deltaT = r/1000;
        
        % This is the listing of how many residuals the different parts of
        % the solution will contribute with.
        % The range gives one residuals per SV and per sample.
        % The velocity gives one residuals per SV and sample (N-1).
        % The acceleration gives one residual per sample (N-2).
        % The jerk gives four residuals (one for the clock) per
        % sample (N-3).
        rNrange = nSv*N;
        rNvel   = nSv*(N-1);
        rNacc   = N-2;
        rNjerk  = 4*(N-3);
        rNwuv   = N-3;      % wind-up-velocity
        
        rNtot = rNrange +rNvel +rNacc + rNjerk + rNwuv;
        
        % Set the weights. We want to rescala them so that the mean range
        % and velocity weight is 1. The difference between them is then
        % completely controlled by the secondary weighting.
        wRanSca = mean(nanmean(wRange));
        wDopSca = mean(nanmean(wDoppler));
        
        w = zeros(rNtot,1);
        for n = 1:nSv
            w(n:nSv:rNrange) = wRange(:,n)'/wRanSca*weights(1);
            w(rNrange+(n:nSv:rNvel)) = wDoppler(:,n)'/wDopSca*weights(2);
        end
        w(rNrange+rNvel+(1:rNacc)) = weights(3);
        
        % Set the weight to zero when the accelerometers might have
        % saturated.
        tempInd = rNrange+rNvel+(1:rNacc);
        w(tempInd(max(abs(accBody),[],2)>15)) = 0;
        
        w(rNrange+rNvel+rNacc+1:end) = weights(4);         % Weighting for the jerk.
        w(rNrange+rNvel+rNacc+1+3*(N-3):end) = weights(5); % Clock error jerk.
        
        w(rNrange+rNvel+rNacc+rNjerk+1:end) = weights(6);  % Wind-up velocity
        
        % Call the 'turnOn' method to only enable some parts of the
        % problem.
        nav_residuals('turnOn', turnOn);
        
        % This is the satellites coordinates in the ECI frame. This is
        % keept to facilitade fixed point iteration. I think the ECI velocity
        % should be used, as opposed to the ECEF.
        svPosCorr = svPos;
        svVelCorr = svVel;
        
        % ##### START OF JACOBIAN #####
        % Start building the index vectors for the sparse Jacobian.
        
        % ##### RANGE #####
        jResTemp = zeros([N,nSv,4]);
        for i = 1:nSv;
            jResTemp(:,i,:) = repmat(N*(i-1)+(1:N)',[1,1,4]);
        end
        
        jParTemp = zeros([N,nSv,4]);
        for i = 1:4;
            jParTemp(:,:,i) = repmat(N*(i-1)+(1:N)',[1,nSv,1]);
        end
        
        jRes = jResTemp(:);
        jPar = jParTemp(:);
        
        % ##### VELOCITY #####
        jResTemp = zeros([N-1,nSv,8+1]);
        for i = 1:nSv;
            jResTemp(:,i,:) = repmat((N-1)*(i-1)+(1:N-1)',[1,1,9]);
        end
        jResTemp = rNrange + jResTemp;
        
        jParTemp = zeros([N-1,nSv,8+1]);
        for i = 1:4;
            for j = 0:1
                jParTemp(:,:,4*j+i) = repmat(N*(i-1)+(1:N-1)'+j,[1,nSv,1]);
            end
        end
        jParTemp(:,:,9) = repmat(4*N+6+(1:N-1)',[1,nSv,1]); %windup
        
        jRes = [jRes; jResTemp(:)];
        jPar = [jPar; jParTemp(:)];
        
        % ##### ACCELERATION #####
        jResTemp = repmat((1:N-2)',[1,3*5]);
        jResTemp = rNrange + rNvel + jResTemp;
        
        jParTemp = zeros([N-2,3*5]);
        for i = 1:3;
            for j = 0:2
                jParTemp(:,3*j+i) = repmat(N*(i-1)+(1:N-2)'+j,[1,1]);
            end
        end
        for i = 1:6;
            jParTemp(:,9+i) = repmat(4*N+i,[N-2,1]);
        end
        
        jRes = [jRes; jResTemp(:)];
        jPar = [jPar; jParTemp(:)];
        
        % ##### JERK #####
        jResTemp = zeros([N-3,4*4]);
        for i = 1:4;
            for j = 0:3
                jResTemp(:,4*j+i) = (N-3)*(i-1)+(1:N-3)';
            end
        end
        jResTemp = rNrange + rNvel + rNacc + jResTemp;
        
        
        jParTemp = zeros([N-3,4*4]);
        for i = 1:4;
            for j = 0:3
                jParTemp(:,4*j+i) = N*(i-1)+(1:N-3)'+j;
            end
        end
        
        jRes = [jRes; jResTemp(:)];
        jPar = [jPar; jParTemp(:)];
        
        % ##### WIND-UP ACC #####
        jResTemp = repmat((1:N-3)',[1,1,3]);
        jResTemp = rNrange + rNvel + rNacc + rNjerk + jResTemp;
        
        jParTemp = zeros([N-3,3]);
        jParTemp(:,1) = 4*N+6+(1:N-3)'; %windup
        jParTemp(:,2) = 4*N+6+(2:N-2)'; %windup
        jParTemp(:,3) = 4*N+6+(3:N-1)'; %windup
        
        jRes = [jRes; jResTemp(:)];
        jPar = [jPar; jParTemp(:)];
        
        % ##### END OF JACOBIAN #####
        
        % Return a simple initial guess.
        beta = zeros(5*N+5,1);
        init = [GET_REC_DATA(0,'pos'),2.5e7];
        for i = 0:3
            beta(N*i+(1:N)) = init(i+1);
        end
        beta(4*N+(1:6)) = [1;1;1; 0;0;0]; %acc slope and bias
        beta(4*N+6+(1:N-1)) = 0; %windup
        
        varargout = {beta, recTime};
    case 'turnOn'
        turnOn  = varargin{1};
        
        wTemp = w;
        wTemp(rNrange+(1:rNvel)) = w(rNrange+(1:rNvel))*(turnOn(1)==1);             % Velocity
        wTemp(rNrange+rNvel+(1:rNacc)) = w(rNrange+rNvel+(1:rNacc))*(turnOn(2)==1); % Acceleration
        wTemp(rNrange+rNvel+rNacc+(1:rNjerk)) = w(rNrange+rNvel+rNacc+(1:rNjerk))*(turnOn(3)==1); % Jerk
        wTemp(end-rNwuv+1:end) = w(end-rNwuv+1:end)*(turnOn(1)==1);                 % Wind-up accelerations
        
        sqrW = sparse(1:rNtot,1:rNtot,sqrt(wTemp));
    case 'Func'
        beta  = varargin{1};
        
        C = 299792458;
        OmegaE = 7.2921151467e-5;
        
        R = zeros(rNtot,1);
        
        % tTransit is the estimated recevier-SV transit time. Used to rotate the
        % SVs back, to use the ECI frame.
        sqrSum = zeros([N,nSv,1]);
        for i = 0:2
            rDiff = repmat(beta(i*N+(1:N)),[1,nSv,1])-svPosCorr(:,:,i+1);
            sqrSum = sqrSum + rDiff.^2;
        end
        tTransit = squeeze(sqrt(sqrSum))/C;
        
        % Rotate the SVs back by the Earths rotation during the transit
        % time, to get the positions in the ECI frame.
        phi = tTransit*OmegaE;
        sinPhi = sin(phi);
        cosPhi = cos(phi);
        
        svPosCorr(:,:,1) =  cosPhi.*svPos(:,:,1) + sinPhi.*svPos(:,:,2);
        svPosCorr(:,:,2) = -sinPhi.*svPos(:,:,1) + cosPhi.*svPos(:,:,2);
        
        phi = (phi(1:end-1,:)+phi(2:end,:))/2;
        sinPhi = sin(phi);
        cosPhi = cos(phi);
        
        svVelCorr(:,:,1) =  cosPhi.*svVel(:,:,1) + sinPhi.*svVel(:,:,2);
        svVelCorr(:,:,2) = -sinPhi.*svVel(:,:,1) + cosPhi.*svVel(:,:,2);
        
        % Read out the accelerometer scale and bias.
        k = beta(4*N+(1:3));
        b = beta(4*N+(4:6));
        
        % ##### START OF JACOBIAN #####
        % Start building the index vectors for the sparse Jacobian.
        
        % ##### RANGE #####
        dt = beta(3*N+(1:N));
        rho = rangeCorr*C;
        
        delPos = zeros([N,nSv,3]);
        for i = 1:3
            delPos(:,:,i) = repmat(beta((i-1)*N+(1:N)),[1,nSv,1])-svPosCorr(:,:,i);
        end
        
        geoRange = squeeze(sqrt(sum(delPos.^2,3)));
        
        res = rho + repmat(dt,[1,nSv]) - geoRange;
        
        R(1:rNrange) = res(:);
        
        % uneVec is the derivative of the geoRange. It is also a unit
        % vector from the receiver to the SV.
        uneVec = -delPos./repmat(geoRange,[1,1,3]); %[N,nSv,3]
        dRdt = ones([N,nSv,1]);
        
        jValTemp = cat(3,uneVec,dRdt); % [N,nSv,4]
        
        jVal = jValTemp(:);
        
        % ##### VELOCITY #####
        d = dopplerCorr; %[N,nSv]
        
        dRec = -diff(dt)/deltaT;
        recVel = diff(reshape(beta(1:3*N),[N,3]))/deltaT;
        windup = repmat(beta(4*N+6+(1:N-1)),[1,nSv]);
        
        % The unit vector from RANGE needs to be resampled.
        u = (uneVec(1:end-1,:,:)+uneVec(2:end,:,:))/2;
        
        relVel = (repmat(reshape(recVel,[N-1,1,3]),[1,nSv,1])-svVelCorr);
        res = d + repmat(dRec,[1,nSv]) + windup - sum(u.*relVel,3);
        
        R(rNrange+(1:rNvel)) = res(:);
        
        dRdX   = +u/deltaT;
        dRdDt  = +ones([N-1,nSv,1])/deltaT;
        
        dRdXp  = -u/deltaT;
        dRdDtp = -ones([N-1,nSv,1])/deltaT;
        
        dRdw   = ones([N-1,nSv,1]);
        
        jValTemp = cat(3,dRdX,dRdDt,dRdXp,dRdDtp,dRdw); % [N-1,nSv,8]
        
        jVal = [jVal; jValTemp(:)];
        
        % ##### ACCELERATION #####
        acc = diff(reshape(beta(1:3*N),[N,3]),2)/deltaT^2;
        
        % Get the Earth's gravity vector in the ECEF frame.
        [gx,gy,gz] = xyz2grav(beta(2:N-1),...
            beta(N+(2:N-1)),...
            beta(2*N+(2:N-1)),...
            deltaT);
        
        % Transform the GPS acceleration to inertial acceleration.
        accIner = acc - [gx,gy,gz];
        
        kRep = repmat(k',[N-2,1]);
        bRep = repmat(b',[N-2,1]);
        
        accBodyMag = repmat(sqrt(sum((accBody.*kRep+bRep+accCorr).^2,2)),[1,3]);
        accInerMag = repmat(sqrt(sum(accIner.^2,2)),[1,3]);
        
        R(rNrange+rNvel+(1:rNacc)) = accBodyMag(:,1)-accInerMag(:,1);
        
        dadX    = -accIner./(accInerMag*deltaT^2);
        dadXP   = 2*accIner./(accInerMag*deltaT^2);
        dadXPP  = -accIner./(accInerMag*deltaT^2);
        
        dadk    = accBody.*(accBody.*kRep+bRep+accCorr)./accBodyMag;
        dadb    = (accBody.*kRep+bRep+accCorr)./accBodyMag;
        
        jValTemp = cat(2,dadX,dadXP,dadXPP,dadk,dadb); % [N-1,nSv,8]
        
        jVal = [jVal; jValTemp(:)];
        
        % ##### JERK #####
        jerk = diff(reshape(beta(1:4*N),[N,4]),3)/deltaT^3;
        
        R(rNrange+rNvel+rNacc+(1:rNjerk)) = jerk(:);
        
        djdX    = -1*ones(N-3,4)/deltaT^3;
        djdXP   = +3*ones(N-3,4)/deltaT^3;
        djdXPP  = -3*ones(N-3,4)/deltaT^3;
        djdXPPP = +1*ones(N-3,4)/deltaT^3;
        
        jValTemp = cat(2,djdX,djdXP,djdXPP,djdXPPP);
        
        jVal = [jVal; jValTemp(:)];
        
        % ##### WIND-UP ACC #####
        angAcc = diff(beta(4*N+6+(1:N-1)),2)/deltaT^2;
        
        R(rNrange+rNvel+rNacc+rNjerk+(1:rNwuv)) = angAcc;
        
        dwvdw  = +1*ones(N-3,1)/deltaT^2;
        dwvdwP = -2*ones(N-3,1)/deltaT^2;
        dwvdwPP = +1*ones(N-3,1)/deltaT^2;
        
        jValTemp = cat(2,dwvdw,dwvdwP,dwvdwPP);
        
        jVal = [jVal; jValTemp(:)];
        
        % ##### BUILD SPARSE JACOBIAN #####
        J = sparse(jRes,jPar,jVal);
        
        % ##### START OF JACOBIAN #####
        
        % NaNs don't contribute to the residuals.
        R(isnan(R)) = 0;
        J(isnan(J)) = 0;
        
        % ### START PLOTTING ###
        if true
            pos = reshape(beta(1:3*N),[N,3]);
            meanPos = nanmean(pos);
            velTime = (recTime(1:end-1)+recTime(2:end))/2;
            s=10;       % a down sampling factor
            
            figure(1)
            
            subplot(2,2,1)  % IMU and GPS acceleration
            accTime = (recTime(1:end-2)+recTime(3:end))/2;
            plot(accTime(1:s:end),accBodyMag(1:s:end,1),'.')
            hold on
            plot(accTime(1:s:end),accInerMag(1:s:end,1),'r')
            ylim([0,30])
            xlabel('Time [sec]')
            ylabel('Acceleration [m/s^2]');
            title('Acceleration magnitude')
            legend('Accelerometers','GPS')
            hold off
            
            subplot(2,2,2)  % Vector valued accelerations
            plot(accTime(1:s:end),accBody(1:s:end,:),'k')
            hold on
            plot(accTime(1:s:end),project_vel(accIner(1:s:end,:),meanPos))
            ylim([-20,20])
            xlabel('Time [sec]')
            ylabel('Acceleration [m/s^2]');
            title('Acceleration vectors')
            legend('Body X', 'Body Y', 'Body Z',...
                'E/W','N/S','Vert')
            hold off
            
            subplot(2,2,3)  % LO offset
            plot(velTime(1:s:end),dRec(1:s:end)*5.2550);
            xlabel('Time [sec]')
            ylabel('LO freq [Hz]');
            title('LO offset')
            
            subplot(2,2,4)  % Wind-up
            plot(velTime(1:s:end),windup(1:s:end,1)*5.2550)
            xlabel('Time [sec]')
            ylabel('Wind-up [Hz]');
            
            %         refAre = 0.45;    % Drag stuff
            %         refAre = pi*(58e-3)^2;
            %         dragConst = 2*1.06/(1.225*refAre);
            %         velAdj = (recVel(1:end-1,:)+recVel(2:end,:))/2;
            %         velMag = sqrt(sum(velAdj.^2,2));
            %         accProj = sum(accIner.*velAdj,2)./velMag;
            %         plot(accTime,-dragConst*accProj./velMag.^2)
            %         xlabel('Time [sec]')
            %         ylabel('Drag coefficient');
            %         title('Drag')
            
            
            figure(2)
            
            subplot(2,2,2) % Vector velocity
            plot(velTime(1:s:end),project_vel(recVel(1:s:end,:),meanPos));
            xlabel('Time [sec]')
            ylabel('Velocity [m/s]');
            title('Velocity local plane')
            legend('E/W','N/S','Vert')
            
            subplot(2,2,4)  % Velocity
            plot(velTime(1:s:end),sqrt(sum(recVel(1:s:end,:).^2,2)));
            xlabel('Time [sec]')
            ylabel('Velocity [m/s]');
            title('Velocity')
            
            subplot(2,2,1)  % Local position
            plot(recTime(1:s:end),project_to_surface(pos(1:s:end,:)));
            xlabel('Time [sec]')
            ylabel('Position [m/s]');
            title('Position local plane')
            
            subplot(2,2,3)  % Position (top-down)
            klm = project_to_surface(pos);
            %         scatter(klm(1:s:end,1),klm(1:s:end,2),[],klm(1:s:end,3),'.')
            plot(klm(1:s:end,1),klm(1:s:end,2))
            axis square
            title('Plane position')
            
            drawnow
        end
        % ### END PLOTTING ###
        
        % Apply weighting to the output.
        if nargout==2
            F = sqrW*R;
            Jw = sqrW*J;
            varargout = {F,Jw};
            
        elseif nargout==0
            CONSTANTS
            
            % Put together som understandable output.
            resHist = reshape(R(1:rNrange),[N,nSv]);
            velResHist = reshape(R(rNrange+(1:rNvel)),[N-1,nSv]);
            
            betaMat = reshape(beta(1:4*N),[N,4]);
            posHist = betaMat(:,1:3);
            velHist = diff(posHist)/deltaT;
            
            Adiag = diag(J'*sqrW.^2*J);
            dopHist = reshape(Adiag(1:4*N),[N,4]);
            
            velDopHist = sqrt(dopHist(1:end-1,:).^2+dopHist(2:end,:).^2)/deltaT;
            
            gpsTimeHist = prRefTime + recTime + beta(3*N+(1:N))/C;
            
            velRecTime = (recTime(1:end-1)+recTime(2:end))/2;
            velDopHist = dopHist;
            
            loHist = -diff(betaMat(:,4))/deltaT*L1/C*2*pi;
            
            accSlope = beta(4*N+(1:3));
            accBias  = beta(4*N+(4:6));
            
            windUp = beta(4*N+6+(1:N-1))*L1/C;
            
            save(sprintf('%sposLS',TRACK_DIRECTORY),'posHist','gpsTimeHist',...
                'recTime','dopHist','resHist', 'velHist','loHist','velRecTime',...
                'velDopHist','velResHist', 'accSlope','accBias','windUp');
        end
end