% Example callback function for analog input logging. Converted from csv to binary using Matlab example logData.m
%   Copyright 2011 The MathWorks, Inc.
% This funcion is called once the number of acquired time points exceeds
%   src.NotifyWhenDataAvailableExceeds.
% The latest analog input data will be written to the log .bin file.
% Add the time stamp and the data values to data. To write data sequentially,
% transpose the matrix.
function logDigData(src, event, fid)
    % Now actually write the results to the log file.
    % Version 1/ms: data = [round(event.TimeStamps*1e3), bi2de(event.Data, 'left-msb')]'; % Convert Digital data into an integer, also store Timestamps as ms. Flip binary TTLs so that channel 0/1 is lowest highest value, consistent with other code orientation
    % Version 2: tenths of ms for sampling rate of 10k
    % Version 2.1: Only save first n-1 channels/ignore final columns as dummy analog input data is last channel: Workaround for disabling analog outputs, to keep timer going: EK 2019/07/26
    data = [round(event.TimeStamps*10e3), bi2de(event.Data(:, 1:end-1), 'left-msb')]'; % Convert Digital data into an integer, also store Timestamps as ms. Flip binary TTLs so that channel 0/1 is lowest highest value, consistent with other code orientation
    fwrite(fid,data,'uint32'); % double is inefficient for digital data but used for time data which is ms as double?
end    
    
% % Original code = saving as .csv: Slow and inefficient
% % Example callback function for analog input logging.
% % This funcion is called once the number of acquired time points exceeds
% %   src.NotifyWhenDataAvailableExceeds.
% % The latest analog input data will be written to the log .csv file.
% % Use plotLogFile() to view the logged data.
% function logAIData(src, event, filename)    
%     % Plot the latest data
%     %plot(event.TimeStamps, event.Data);
%     
%     % Now actually write the results to the .csv log file.
%     dlmwrite(filename,[ event.TimeStamps event.Data],'delimiter',',','-append');
