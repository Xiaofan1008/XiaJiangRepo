function loopedProcessIntan
base = 'X:\Tim';
%% List of directories to process
listOfDirectories = [
    %'\RAT0007\';...
    %'\RAT0009\';...
    '\RAT0010\';... % Up to 18 of 19
    %'\RAT0011\';...
    %'\RAT0015\';...
    %'\RAT0016\';...
    %'\RAT0017\';...
    ];
%% Loop through each directory
tic;
for d = 1:1:size(listOfDirectories,1)
    cd([base listOfDirectories(d,:)]);
    % Loop through each subfolder
    files = dir;
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    for f = 18:length(subFolders)        
        fprintf('Now processing . . . %d of %d\n',(f-2),length(subFolders)-2);
        cd([base listOfDirectories(d,:) subFolders(f).name]);
        dnName = dir('*dn.dat');
        if ~isempty(dnName)
            delete(dnName.name);
        end
        dnName = dir('*mu.dat');
        if ~isempty(dnName)
            delete(dnName.name);
        end
        try
            processIntan(1,0);
        catch
            disp('Error. Bad filename');
        end
    end
    fprintf('\n');
    fprintf('Loop %d of %d complete!\n',d,size(listOfDirectories,1));
    toc;
    fprintf('\n');
end
end
