function processIntan_h
dbstop on error
% Step One: Remove CAR and artefact
ifile = pwd;
denoiseIntan(ifile);
% Step Two: Generate mu and lfp files
allExtract;
% Step Three: Generate spikes
splitSpikes;
end