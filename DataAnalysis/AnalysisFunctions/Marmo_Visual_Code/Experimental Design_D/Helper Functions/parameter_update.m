%% Script for sending a parameter update over TCP/IP

function parameter_update(t,params,n)
% Convert the parameter set into strings
CHANNEL = strcat('channel=',(params{1}));
pT = strcat('pulseOrTrain=',num2str(params{7}));
nPT = strcat('numberOfStimPulses=',num2str(params{8}));
pTP = strcat('pulseTrainPeriod=',num2str(params{9}));
fPD = strcat('firstPhaseDuration=',num2str(params{13}));
sPD = strcat('secondPhaseDuration=',num2str(params{14}));
IPD = strcat('interphaseDelay=',num2str(params{15}));
fPA = strcat('firstPhaseAmplitude=',num2str(params{16}));
sPA = strcat('secondPhaseAmplitude=',num2str(params{17}));
% prSD = strcat('postTriggerDelay=',num2str(params{19}+10));
prAS = strcat('preStimAmpSettle=',num2str(params{19}));
prSD = strcat('postTriggerDelay=',num2str(params{6}+params{19}+10));
poSD = strcat('refractoryPeriod=',num2str(params{20}+100));
poAS = strcat('postStimAmpSettle=',num2str(params{20}));
prCR = strcat('postStimChargeRecovOn=',num2str(params{23}));
poCR = strcat('postStimChargeRecovOff=',num2str(params{24}));
polarity = strcat('stimPolarity=',num2str(params{12})); %cathodic=0, anodic=1

% Send the parameters over the TCP port
fprintf(t,CHANNEL);
if strcmp(fPA,'firstPhaseAmplitude=-1') || strcmp(fPD,'firstPhaseDuration=-1')
    fprintf(t,'enabled=0');
else
    fprintf(t,'enabled=1');
end
fprintf(t,pT);
fprintf(t,nPT);
fprintf(t,pTP);
fprintf(t,fPD);
fprintf(t,sPD);
fprintf(t,IPD);
fprintf(t,fPA);
fprintf(t,sPA);
fprintf(t,'stimShape=1'); %'Biphasic'=0, 'BiphasicWithInterphaseDelay'=1, 'Triphasic'=2
fprintf(t,polarity);
fprintf(t,prSD);
fprintf(t,prAS);
fprintf(t,poSD);
fprintf(t,poAS);
fprintf(t,'enableAmpSettle=1');
fprintf(t,'maintainAmpSettle=1');
fprintf(t,prCR);
fprintf(t,poCR);
fprintf(t,'enableChargeRecovery=1');
fprintf(t,['ID=' num2str(n)]);
end


