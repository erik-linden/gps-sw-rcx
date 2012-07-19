% Decodes FFU Main FPGA memory.
% Input data must be in binary format, alternating data and time packets.
% The file must start with a data packet for MUX channel 0. (Use hex editor)
% The file should end with a data packet for channel 0xf followed by its
% time packet.


function [adc1, adc2, status, time, length] = read_decoded(filename)
  fid = fopen(filename);
  
  clear adc1;
  clear adc2;
  clear time;
  clear length;
    
  signature = fread(fid, 4, '*uint8', 0, 'b')';
  if signature ~= 'FFD0'
       disp(['    Error: Missing signature.']);
       fclose(fid);
       return;
  end
  
  for muxch = 1:16
    datapoints = fread(fid, 1, 'uint32', 0, 'b');
        
    A = fread(fid, [3 datapoints], '*uint32', 0, 'b');
    
    time(muxch, 1:datapoints) = double(A(1, :)) / 1000;
    status(muxch, 1:datapoints) = A(2, :);
    
    adcdata = A(3, 1:datapoints);
    adc1(muxch, 1:datapoints) = single(bitand(adcdata, hex2dec('ffff0000')) / hex2dec('10000')); % mask and shift
    adc2(muxch, 1:datapoints) = single(bitand(adcdata, hex2dec('0000ffff'))); % mask
    
    length(muxch, 1) = datapoints;
  end

  fclose(fid);
end
