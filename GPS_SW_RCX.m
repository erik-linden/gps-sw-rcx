function GPS_SW_RCX()

clear global
CONSTANTS_H;
CONSTANTS;

choose = 0;
while(choose>4 || choose<1)
    fprintf('Would you like to:\n1) Track 1 or more satellites\n2) Obtain LS navigation solution\n3) Obtain legacy navigation solution\n4) Quit\n')
    choose = input(': ');
end

switch choose
    case  1
        prnFlag = 1;
        while(prnFlag)
            prn = input('\nPlease enter one or more satellites to track in the form [SV1 SV2 ...]\n or press enter to search for all satellites: ');
            if(isempty(prn))
                prn = 1:32;
                prnFlag =0;
            else
                prnFlag = 0;
            end
        end
        
        %the number of seconds to track the signal for
        Nfiles = input('\nEnter the number of seconds to track the signal for: ');
        
        fileNo = input('\nEnter which file to start with: ');fprintf('\n');
        
        tStart = tic;
        for x=1:length(prn)
            % Start by loading the first file.
            in_sig = LOAD_GPS_DATA(FILE_NAME,fileNo);
            
            % Set the time offset of the file. Assume all files have the same
            % length.
            t0 = fileNo*length(in_sig)*TP;
            
            % Figure out where to place the center frequency. This is used for
            % acquisition.
            if USE_AIDING
                [df, elev] = ESTIMATE_DOPPLER(t0, prn(x));
                dLO = GET_REC_DATA(t0, 'dFre');
                F0 = df - dLO;
                string = sprintf(' Elevation: %+4.f, ',elev);
            else
                F0 = 0;
                elev = [];
                string = '';
            end
            
            % Do aquisition. The coherent integration time is fixed to 20
            % milliseconds.
            if ~USE_AIDING || elev>-20
                code = CACODEGN(prn(x));
                corr_aq=DIGITIZE_CA(code,0,N_CODES_AQU*ONE_MSEC_SAM);
                [df, cst_k, magnitude] = FFT_ACQUISITION(in_sig, corr_aq, F0);
                
                % If the signal was not found, quit this satellite.
                if(CNO(magnitude,N_CODES_AQU*1e-3) < DB_DETECTION)
                    fprintf('PRN %02d NOT FOUND%s CNo: %2.2f\n', prn(x), string, CNO(magnitude,N_CODES_AQU*1e-3));
                    % Otherwise, track the satellite.
                else
                    fprintf('PRN %02d FOUND    %s CNo: %2.2f, Doppler Frequency: %+3.f\n', prn(x), string, CNO(magnitude,N_CODES_AQU*1e-3), df-F0);
                    if Nfiles>0
                        corr = GEN_CORR(code);
                        SIGNAL_TRACKING(df, cst_k, in_sig, prn(x), corr, FILE_NAME, fileNo, t0, Nfiles)
                    end
                end
            else
                fprintf('PRN %02d IGNORED  %s\n', prn(x), string(1:end-2));
            end
        end
        fprintf('Elapsed time is %s\n',datestr(datenum(0,0,0,0,0,toc(tStart)),'HH:MM:SS'))
    case 2
        sv = LIST_AVAIL_TRACKS();
        svInp = input(['Enter satellites to obtain navigation solution from or \n'...
            'press enter for [' sprintf('%02u ',sv) ']\n:']);
        if(~isempty(svInp))
            sv = svInp;
        end
        
        disp('')
        disp('Enter solver controls as a N-by-4 matrix with the following format:')
        disp('|R1 V1 A1 J1|')
        disp('|R2 V2 A2 J2|')
        disp('| .  .  .  .|')
        disp('|RN VN AN JN|')
        disp('Each row is one iteration of the solver. R is a downsampling factor.')
        disp('V, A and J are booleans which enables or disables Doppler, accelerometers')
        disp('and jerk, respectively.')
        disp('Press enter to use default:')
        disp('[1000 0 0 1;...')
        disp(' 100  1 0 1;...')
        disp(' 10   1 0 1];')
        
        setInp = input(':');
        if(isempty(setInp))
            setInp = [1000 0 0 1;...
                      100  1 0 1;...
                      10   1 0 1];
        end
        
        GLOBAL_LS_SOL(sv,setInp)
    case 3
        sv = LIST_AVAIL_TRACKS();
        svInp = input(['Enter satellites to obtain navigation solution from or \n'...
            'press enter for [' sprintf('%02u ',sv) ']\n:']);
        if(~isempty(svInp))
            sv = svInp;
        end
        
        rPos = input('Please enter position sampling reduction factor (1 => 1000 Hz): ');
        useIcp = input('Use differential mode? (default no): ');
        rVel = input('Please enter velocity sampling reduction factor (1 => 1000 Hz): ');
        
        % Solve for position.
        [posHist, gpsTimeHist, recTime, dopHist, resHist]=POSITION_TRACK(sv, rPos, useIcp);
        
        % Solve for velocity.
        [velHist loHist velRecTime velDopHist velResHist] = VELOCITY_TRACK(sv, posHist, recTime, gpsTimeHist, rVel);
        
        % Save data.
        save('pos','posHist','gpsTimeHist','recTime','dopHist','resHist',...
            'velHist','loHist','velRecTime','velDopHist','velResHist');
    otherwise
        return
end