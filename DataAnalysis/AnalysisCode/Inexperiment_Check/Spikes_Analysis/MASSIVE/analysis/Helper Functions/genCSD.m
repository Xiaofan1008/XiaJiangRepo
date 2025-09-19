function genCSD(path)
% This function uses a flash datafile to set up a CSD matrix and preps CSDplotter
if nargin
    cd(path);
else
    path = pwd;
end
%% Variables
BN = [-100, 400]; % Window in ms
DEAD = [2,9]; % Dead channels
%% Loading
try
    lfp = loadLFP;
catch
    path = uigetdir;
    cd(path)
    lfp = loadLFP;
end
trig = loadTrig(1);
DEPTH = Depth;
nTrig = length(trig);
nChn = size(lfp,1);
nDead = length(DEAD);
LFP = zeros(nChn,diff(BN)+1);
%% Logic
for t = 1:nTrig
    for c = 1:nChn
        LFP(c,:) = LFP(c,:) + lfp(DEPTH(c),trig(t)+BN(1):trig(t)+BN(2));
    end
end
LFP = LFP ./ nTrig;
%% Remove Dead Channels
for d = 1:nDead
    LFP(DEAD(d),1) = 9e5;
end
LFP(LFP(:,1) == 9e5,:) = [];
figure; hold on;
%% Saving
t_CSDname = strsplit(path,'\');
t_CSDname = [t_CSDname{end} '.mat'];
t_CSDpath = 'C:\Users\tall0003\Google Drive\Monash PhD\Scripts\Analysis Tools\toolboxes\csd\CSDplotter-0.1.1\Data\';
% Invert the CSD file for CSDPlotter
LFP = flipud(LFP);
save([t_CSDpath t_CSDname],'LFP');
LFP = flipud(LFP);
%% Run a Spline CSD
dt = 1; % Resolution of the LFP in ms
runCSD(LFP,dt,32,DEAD);
YLIM = ylim;
line([abs(BN(1)) abs(BN(1))],[YLIM(1) YLIM(2)],'Color','k','LineWidth',2);
%text(abs(BN(1))+1,YLIM(2)-50,'Stimulus','Color','k','FontSize',20);
%% Closing
end