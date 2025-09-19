folder = dir;
names = {folder([folder.isdir]).name};
here = pwd;
for i = 3:length(names)
    path = names{i};
    file = dir(path);
    subnames = {file([file.isdir]).name};
    for j = 3:length(subnames)
        filepath = subnames{j};
        finalpath = [here '\' path '\' filepath];
        documentRecovery(finalpath);
        disp([finalpath ': DONE.']);
    end
end
close all
clearvars
