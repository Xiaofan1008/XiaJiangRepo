%% Function uses a high-pass filter to detrend
function filteredwave = Detrend_Short(wave,HZ,FS,ORDER)
% Create a filter form the provided parameters
[b,a] = butter(ORDER,HZ/(FS/2),'high');
filteredwave = filtfilt(b,a,wave);
end