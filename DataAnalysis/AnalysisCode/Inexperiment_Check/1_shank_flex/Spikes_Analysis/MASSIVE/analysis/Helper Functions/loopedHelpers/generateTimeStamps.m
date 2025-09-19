%% Script that takes a directory and generates timestamps from the digitalin
fileinfo = dir([filepath, '\digitalin.dat']);
if ~isempty(fileinfo)
    num_samples = fileinfo.bytes/2;
    digin_fid = fopen([filepath, '\digitalin.dat'],'r');
    % I should note that the digitalin.dat file in this case has two rows,
    % one for electrical stimulation, and one for visual stimulation
    digital_in = fread(digin_fid, num_samples, 'uint16');
    fclose(digin_fid);
    stimDig = flip(find(digital_in == 1)); % Fix for finding the falling line instead of the rising line
    visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
    dt = diff(stimDig);
    kill = dt == -1;
    stimDig(kill) = [];
    dt = diff(visDig);
    kill = dt == -1;
    visDig(kill) = [];
    nStamps = max([length(stimDig); length(visDig)]);
    time_stamps = nan(2,nStamps);
    time_stamps(1,1:length(stimDig)) = flip(stimDig);
    time_stamps(2,1:length(visDig)) = flip(visDig);
    time_stamps = time_stamps ./ (FS/1e3); % Convert to milliseconds
    if isnan(time_stamps(2,2))
        time_stamps = time_stamps(1,:);
    else
        time_stamps = time_stamps(2,:);
    end
end

