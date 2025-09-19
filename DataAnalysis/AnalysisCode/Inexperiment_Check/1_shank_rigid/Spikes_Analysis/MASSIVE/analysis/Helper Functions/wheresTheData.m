%% Script for computing which recordings have actual data in them
path = 'E:\CJ196 Intan';
RECORDINGS = dir(path);
RECORDINGS(1:2) = [];
nDATA = zeros(1,length(RECORDINGS));
for i = 1:length(RECORDINGS)
    datafile = dir([path '\' RECORDINGS(i).name '\*_exp_datafile*.mat']);
    if isempty(datafile)
        nDATA(i) = 0;
    else
        nDATA(i) = 1;
    end
end
RECORDINGS(nDATA == 0) = [];
%% Example code for indexing through recordings
% for i = 4:length(RECORDINGS)
% allSpikes([path '\' RECORDINGS(i).name]);
% end

path = 'E:\CJ197 Intan';
RECORDINGS = dir(path);
RECORDINGS(1:2) = [];
nDATA = zeros(1,length(RECORDINGS));
for i = 1:length(RECORDINGS)
    datafile = dir([path '\' RECORDINGS(i).name '\*_exp_datafile*.mat']);
    if isempty(datafile)
        nDATA(i) = 0;
    else
        nDATA(i) = 1;
    end
end
RECORDINGS(nDATA == 0) = [];
%% Example code for indexing through recordings
for i = 1:length(RECORDINGS)
allSpikes([path '\' RECORDINGS(i).name]);
end