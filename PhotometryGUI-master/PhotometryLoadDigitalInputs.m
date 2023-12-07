function [data_bin, ts_dig, dig_Fs, dig_file_ver, dig_num]= PhotometryLoadDigitalInputs(filename)


%% Load file
% filename = 'Test_20190719.dig';

fid = fopen(filename, 'rb');
file_data = fread(fid, 'uint32');
fclose(fid);

%% Process data
dig_file_ver = file_data(1);
if dig_file_ver == 1
    dig_num = file_data(2);
    ts_dig = file_data(3:2:end);
    ts_dig = ts_dig / 1e3; % Timestamps saved as ms to save as int, convert to sec
    data_dec = file_data(4:2:end);
    data_bin = de2bi(data_dec, dig_num, 'left-msb');
    dig_Fs = 1 / median(diff(ts_dig)); % Saved as some other unit rather than sec, but converted to sec prev
elseif dig_file_ver == 2
    dig_Fs = file_data(2); % Resolution of timestamps as recorded
    dig_res_us = file_data(3); % Resolution of timestamps as recorded
    dig_num = file_data(4);
    ts_dig = file_data(5:2:end);
    ts_dig = ts_dig / 1e6 * dig_res_us; % Timestamps saved as multiple of us to save as int, convert to sec
    data_dec = file_data(6:2:end);
    if ~isempty(data_dec)
        data_bin = de2bi(data_dec, dig_num, 'left-msb');
    else
        data_bin = [];
    end
else
    fprintf('Unrecognized digital file version = %d ?\n', dig_file_ver);
end


% %% Plot data
% plot(ts_dig, data_dec);

