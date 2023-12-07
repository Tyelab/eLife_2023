function fid = FilePrepDigData(logDigFile, Fs, handles)

fid = fopen(logDigFile, 'wb');
% % Version 1:
% fwrite(fid, 1, 'uint32'); % File version number
% fwrite(fid, sum(strcmp('InputOnly', {handles.s.Channels.MeasurementType})), 'uint32'); % Num digital channels

% Version 2:
fwrite(fid, 2, 'uint32'); % File version number
fwrite(fid, Fs, 'uint32'); % File version number
fwrite(fid, 100, 'uint32'); % Duration of each bin in us/time resolution: 100us = 0.1ms = 10e3 Hz: flexible in case want different resolution than sampling frequency (likely overly complicated, however)
fwrite(fid, sum(strcmp('InputOnly', {handles.s.Channels.MeasurementType})), 'uint32'); % Num digital channels
