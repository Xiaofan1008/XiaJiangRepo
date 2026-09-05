%% Master function for Post-Experiment looped scripts
function master(filepath)
%% User-defined variables
threshsign = 1;
ARTIFACT_CUTOFF = 0.3e3;
BIN = [-200 100];
spWINDOW = [1 11];
% IMPORTANT. Must set the high-level save directory first. Requires no '\'
% on the end.
save_directory = filepath(1:end-35);
%% Basic Analysis
loopedSpikes(filepath,threshsign,ARTIFACT_CUTOFF,BIN,spWINDOW,save_directory);
loopedPSTH(filepath,threshsign,ARTIFACT_CUTOFF,BIN,spWINDOW,save_directory);
loopedResponseCurve(filepath,threshsign,ARTIFACT_CUTOFF,BIN,spWINDOW,save_directory);
end

