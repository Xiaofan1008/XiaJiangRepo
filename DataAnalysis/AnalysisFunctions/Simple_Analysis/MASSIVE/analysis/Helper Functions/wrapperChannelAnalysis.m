function wrapperChannelAnalysis()
% This function is a wrapper for scResponse() in order to save relevant
% plots for each channel
%% Variables
base = 'X:\Tim';
saveBase = 'C:\Users\tall0003\Google Drive\Monash PhD\Data\RAT0017 Advanced\Channel Characterisation';
listOfDirectories = {'\RAT0017\estim_pen1_009_190401_193112';...
    '\RAT0017\estim_pen1_010_190401_201325';...
    '\RAT0017\estim_pen1_011_190401_204831';...
    };
saveDir = {'\Rec_009';...
    '\Rec_010';...
    '\Rec_011';...
    };
nChn = 32;
%% Conversion
REC = listOfDirectories(:,1)';
saveREC = saveDir(:,1)';
nRec = size(REC,2);
%% Logic
for n = 2:nRec
    cd([base REC{n}]);
    for c = 1:nChn
        scResponse(c,0);
        % Save
        saveas(gcf,[saveBase saveREC{n} '\Chn_0' num2str(c) '.png']);
        close all;
    end
end

%% Closing