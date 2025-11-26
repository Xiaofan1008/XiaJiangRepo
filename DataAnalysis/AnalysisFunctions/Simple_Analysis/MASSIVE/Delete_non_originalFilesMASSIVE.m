%% delete all non-original files on massive
folder = fileparts(which('MASSIVE_ANALYSIS_SM')); % Determines filepath to folder containing your .m file.
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'dual_vary_noshift_pen4_002_200707_205943'])%change working directory to where data is stored - currently manually input
%%
D = dir; % A is a struct ... first elements are '.' and '..' used for navigation.
for k = 3:length(D) % avoid using the first ones
    currD = D(k).name; % Get the current subdirectory name
    try
        cd([D(k).folder filesep currD]) % change the directory (then cd('..') to get back)
    catch
        continue
    end
    
    if ~isempty(dir('*_dn_sab.dat'))
        delete *_dn_sab.dat
    end
    
    if ~isempty(dir('*.mu_sab.dat'))
        delete *.mu_sab.dat
    end
    
    if ~isempty(dir('*.sp.mat'))
        delete *.sp.mat
    end
    
    if ~isempty(dir('*trig.dat'))
        delete *trig.dat
    end
    
    if ~isempty(dir('IDstruct.mat'))
        delete IDstruct.mat
    end
    
    if ~isempty(dir('Averagetrialresponse.mat'))
        delete Averagetrialresponse.mat
    end
end
%%
if ~isempty(dir('*_dn_sab.dat'))
    delete *_dn_sab.dat
end

if ~isempty(dir('*.mu_sab.dat'))
    delete *.mu_sab.dat
end

if ~isempty(dir('*.sp.mat'))
    delete *.sp.mat
end

if ~isempty(dir('*trig.dat'))
    delete *trig.dat
end

if ~isempty(dir('IDstruct.mat'))
    delete IDstruct.mat
end

if ~isempty(dir('Averagetrialresponse.mat'))
    delete Averagetrialresponse.mat
end
%%
fprintf('Deleted all non-original files')