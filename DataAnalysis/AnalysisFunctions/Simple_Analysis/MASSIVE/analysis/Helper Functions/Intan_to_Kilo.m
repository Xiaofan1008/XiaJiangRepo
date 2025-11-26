function Intan_to_Kilo(filepath,kfilepath)
%% This function converts an Intan amplifier file into
% a KiloSort compatible format
trig = [];

%% Initial Setup
nChn = 32;
FS = 30000;
SCALEFACTOR = 10;
loadTrig;
trig = cast(trig,'double');
blankingBin = [-100 100];
lSize = 60*(1e3);
iFile = [filepath '\amplifier.dat'];

%% Common Artefact Rejection
applyCARtoDat(iFile,nChn);

%% Setup for loop
iFile = [filepath '\amplifier_CAR.dat'];
iFID = fopen(iFile,'r');
iData = dir(iFile);
iBytes = iData.bytes/(nChn * 2);
kFile = [kfilepath '\kilo_amplifier.dat'];
kFID = fopen(kFile,'w');
nLoops = ceil(iBytes/FS) / (lSize/1e3);

%% Main Logic
for n = 1:nLoops
    tic;
    disp(['Now running loop ' num2str(n) ' of ' num2str(nLoops)]);
    theseTrig = trig(trig > (n-1)*lSize - blankingBin(1)*30 & trig <= (n)*lSize - blankingBin(2)*30);
    iRaw = fread(iFID,[nChn,(lSize * (FS/1e3))],'int16');
    iRaw = iRaw .* 0.195;
    for c = 1:nChn
        for t = 1:length(theseTrig)
            thisTrig = cast(theseTrig(t),'int32');
            thisTrig = thisTrig - ((n-1)*(lSize * (FS/1e3))/30);
            iRaw(c,(thisTrig + blankingBin(1))*30:(thisTrig + blankingBin(2))*30) = BlankArte(iRaw(c,(thisTrig + blankingBin(1))*30:(thisTrig + blankingBin(2))*30),blankingBin);
        end
    end
    iRaw = artefactBlank(iRaw);
    kRaw = iRaw .* SCALEFACTOR;
    fwrite(kFID,kRaw,'int16');
    toc;
end

%% Close Down
fclose(iFID);
fclose(kFID);
end