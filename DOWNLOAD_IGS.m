function [C, ephType, I, ionType] = DOWNLOAD_IGS(week,day,year,dayOfYear)
% [C, ephType, I, ionType] = DOWNLOAD_IGS(week,day,year,dayOfYear)
% 
% Downloads latest eph. and ion files from the IGS service.
% 
% Output:
% 	C	  	-	Eph. structure
% 	ephType	-	Type of eph. (0->None,1->Final,2->Rapid,3->Ultrarapid)
% 	I	  	-	Eph. structure
% 	ionType	-	Type of ion  (0->None,1->Final,2->Rapid)

adress = 'cddis.gsfc.nasa.gov';
ephDir = sprintf('~/gps/products/%04d',week);
ionDir = sprintf('~/gps/products/ionex/%04d/%03d',year,dayOfYear);
dirOpen = 1;
ephType = 0;
ionType = 0;
C = [];
I = [];

fprintf('Searching for ephemerides file.\n')
try
    f = ftp(adress);
    cd(f,ephDir);
catch ME1
    warning(ME1.identifier,'Error while opening ''ftp://%s%s.''',adress,ephDir(2:end));
    dirOpen = 0;
end

if dirOpen
    list=dir(f);
    l = length(list);
    found = 0;
    
    k = 1;
    while k<=2 && ~found
        switch k
            case 1
                % Final igswwwwd.sp3
                nameStr = sprintf('igs%04d%1d.sp3',week,day);
                strLength = 12;
            case 2
                % Rapid igrwwwwd.sp3
                nameStr = sprintf('igr%04d%1d.sp3',week,day);
                strLength = 12;
        end
        n = 1;
        while n<=l && ~found
            if strncmp(nameStr,list(n).name,strLength)
                ephType = k;
                fileN = n;
                found = 1;
            end
            n=n+1;
        end
        k=k+1;
    end
    
    if ~found
        %Ultrarapid iguwwwwd_hh.sp3
        nameStr = sprintf('igu%04d%1d',week,day);
        strLength = 8;
        indList = [];
        
        n = 1;
        while n<=l
            name = list(n).name;
            if strncmp(nameStr,name,strLength) && strcmp('.sp3',name(12:15));
                ephType = 3;
                indList = [indList, n];
                found = 1;
            end
            n=n+1;
        end
        
        hourCurr = -10;
        for n = 1:length(indList)
            k = indList(n);
            strHour = list(k).name(10:11);
            hour = str2double(strHour);
            if hour>=hourCurr
                fileN = k;
            end
        end
    end
    
    if found
        filename = list(fileN).name;
        mget(f,filename);
        
        switch ephType
            case 1
                disp('Found: Final.')
            case 2
                disp('Found: Rapid.')
            case 3
                disp('Found: Ultra-Rapid.')
        end
        
        [notUnp,result] = dos(sprintf('gunzip %s',filename));
        if notUnp
            warning('RAIN:NavData','%s',result(1:end-1))
        end
        C = GEN_EPH_LUT({filename(1:end-2)});
        delete(filename(1:end-2))
        disp('Ephemereis LUT created.')
    else
        warning('RAIN:NavData','No ephemereis file found.')
    end
end


% Ion
fprintf('\nSearching for ionex file.\n')

dirOpen = 1;
try
    f = ftp(adress);
    cd(f,ionDir);
catch ME1
    warning(ME1.identifier,sprintf('Error while opening ftp://%s%s.',adress,ionDir(2:end)));
    dirOpen = 0;
end

if dirOpen
    list=dir(f);
    l = length(list);
    found = 0;
    yearShort = year-floor(year/100)*100;
    
    k = 1;
    while k<=2 && ~found
        switch k
            case 1
                % Final igsgddd0.yyi
                nameStr = sprintf('igsg%03d0.%02di',dayOfYear,yearShort);
                strLength = 12;
            case 2
                % Rapid igrgddd0.yyi
                nameStr = sprintf('igrg%03d0.%02di',dayOfYear,yearShort);
                strLength = 12;
        end
        n = 1;
        while n<=l && ~found
            if strncmp(nameStr,list(n).name,strLength)
                ionType = k;
                fileN = n;
                found = 1;
            end
            n=n+1;
        end
        k=k+1;
    end
    
    if found
        filename = list(fileN).name;
        mget(f,filename);
        
        switch ephType
            case 1
                disp('Found: Final.')
            case 2
                disp('Found: Rapid.')
        end
        
        [notUnp,result] = dos(sprintf('gunzip %s',filename));
        if notUnp
            warning('RAIN:NavData','%s',result(1:end-1))
        end
        I = GEN_ION_LUT(filename(1:end-2));
        delete(filename(1:end-2))
        disp('Ion LUT created.')
    else
        fprintf('No ionex file found.\n')
    end
end