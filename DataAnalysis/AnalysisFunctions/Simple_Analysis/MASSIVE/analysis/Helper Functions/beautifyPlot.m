function ax = beautifyPlot(fS,ax)
%% This function takes an axis handle and applies standard parameters to it
if nargin<2
    ax = gca;
end
if ~nargin
    fS = 30;
end
%% Standard Settings
fontSize = fS;
lineWidth = 2;
%% Apply Settings
box off
grid off
set(ax,'FontSize',fontSize,...
       'LineWidth',lineWidth,...
       'TickDir','out');       
end