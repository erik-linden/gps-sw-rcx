function in_sig = LOAD_GPS_DATA(file, fileNo)
% in_sig = LOAD_GPS_DATA(file, fileNo)
% 
% INPUTS
% file      GPS data file, either a *.bin or a parsed *.mat
% fid       The file pointer
% fileNo    The current file number

% OUTPUTS
% in_sig    1 second of +/-1's and +/-3's data
% fid       the file pointer
% fileNo    the next file in the sequence of files

CONSTANTS_H;

if IS_COMPLEX
    file = sprintf('%sI%d.dat',file,fileNo);
    fid = fopen(file);
    if(fid<0)
        fprintf('\nCould not read in-phase file %s',file)
    end
    in_sig_I = fread(fid)-3;
    fclose(fid);
    
    file = sprintf('%sQ%d.dat',file,fileNo);
    fid = fopen(file);
    if(fid<0)
        fprintf('\nCould not read quadrature file %s',file)
    end
    in_sig_Q = fread(fid)-3;
    fclose(fid);
    
    in_sig = in_sig_I+1i*in_sig_Q;
else
    file = sprintf('%s%d.dat',file,fileNo);
    fid = fopen(file);
    if(fid<0)
        fprintf('\nCould not read file %s',file)
    end
    
    in_sig = fread(fid)-3;
    in_sig = in_sig;
    
    fclose(fid);
end

return;