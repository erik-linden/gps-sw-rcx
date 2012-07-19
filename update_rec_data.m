% Replaces the content of the rec.mat file with data generated in the
% current position solution.

CONSTANTS
load(RECEIVER_FILE)

% Read time from rec file.
week = Rec.GPSWeek;
startTime = Rec.GPSTime(1);

if Rec.EphType~=1 || Rec.IonType~=1
    fprintf('\nTrying to update data files.\n')
    % Derive date information.
    ymdhms = gps_to_ymdhms([week, startTime]);
    year = ymdhms(1);
    day = floor(startTime/(3600*24));
    dayOfYear = floor(datenum(ymdhms)-datenum([year 0 0 0 0 0]));
    
    % Download latest data files.
    fprintf('\n')
    [C, ephType, I, ionType] = DOWNLOAD_IGS(week,day,year,dayOfYear);
    
    % Check if data is newer.
    if ephType<Rec.EphType || Rec.EphType==0
        fprintf('\nUpdating eph. from type %d to type %d.',Rec.EphType,ephType)
        Rec.Eph = C;
        Rec.EphType = ephType;
    else
        fprintf('\nEph has not been updated.\n')
    end
    
    if ionType<Rec.IonType || Rec.IonType==0
        fprintf('Updating ion from type %d to type %d.',Rec.IonType,ionType)
        Rec.Ion = I;
        Rec.IonType = ionType;
    else
        fprintf('Ion has not been updated.\n')
    end
else
    fprintf('\nData files are final versions.\n')
end

% See if there is any position data.
if exist('pos.mat','file')
    update = input(['\nA position solution file was found.\n',...
        'Use it to update the receiver file?: ']);
    
    if update
        load('pos')
        
        % Warn if there might be a file missmatch.
        posGpsTime = gpsTimeHist(1);
        if abs(posGpsTime-startTime)>60*10
            errorOk = input(sprintf(['\nThe difference in time between the the\n',...
                'recevier file start time, and the position-\n',...
                'file start time is %.0f minutes.\n',...
                'Continue?: '],...
                (posGpsTime-startTime)/60));
        else
            errorOk = 1;
        end
        
        if errorOk
            fprintf('\nUpdating with position file.\n')
            
            % See what position data is good.
            indGoodP = ~isnan(gpsTimeHist);
            
            % Update the rec structure.
            Rec.Time = recTime(indGoodP)';
            Rec.GPSTime = gpsTimeHist(indGoodP)';
            Rec.Pos = posHist(indGoodP,:)';
            
            % See what velocity data is good.
            indGoodV = ~isnan(loHist);
            
            % Interpolat the velocity.
            velInt = zeros(sum(indGoodP),3);
            for n = 1:3
                velInt(:,n) = interp1(velRecTime(indGoodV),velHist(indGoodV,n),recTime(indGoodP),'linear','extrap');
            end
            
            % Update the velocity data.
            Rec.Vel = velInt';
            
            % Interpolate and update the lo data.
            loInt = interp1(velRecTime(indGoodV),loHist(indGoodV),recTime(indGoodP),'linear','extrap');
            Rec.dFre = loInt'/(2*pi);
        else
            fprintf('\nSkipping position file.\n')
        end
    else
        fprintf('\nIgnoring position file.\n')
    end
end

% Save the new file.
save(RECEIVER_FILE,'Rec')