%% Script for sending a parameter update over TCP/IP
function parameter_update(t,params,n)
% Convert the parameter set into strings
CHANNEL = strcat('channel=',(params{1}));
trEL = strcat('triggerEdgeOrLevel=',num2str(params{4}));
pTD = strcat('postTriggerDelay=',num2str(params{6}));
pT = strcat('pulseOrTrain=',num2str(params{7}));
nPT = strcat('numberOfStimPulses=',num2str(params{8}));
pTP = strcat('pulseTrainPeriod=',num2str(params{9}));
pSR = strcat('postStimRefractory=',num2str(params{10}));
stS = strcat('stimShape=',num2str(params{11}));
stP = strcat('stimPolarity=',num2str(params{12}));
fPD = strcat('firstPhaseDuration=',num2str(params{13}));
sPD = strcat('secondPhaseDuration=',num2str(params{14}));
IPD = strcat('interphaseDelay=',num2str(params{15}));
fPA = strcat('firstPhaseAmplitude=',num2str(params{16}));
sPA = strcat('secondPhaseAmplitude=',num2str(params{17}));
enAS = strcat('enableAmpSettle=',num2str(params{18}));
prAS = strcat('preStimAmpSettle=',num2str(params{19}));
poAS = strcat('postStimAmpSettle=',num2str(params{20}));
maAS = strcat('maintainAmpSettle=',num2str(params{21}));
enCR = strcat('enableChargeRecovery=',num2str(params{22}));
prCR = strcat('postStimChargeRecovOn=',num2str(params{23}));
poCR = strcat('postStimChargeRecovOff=',num2str(params{24}));

% Send the parameters over the TCP port
fprintf(t,CHANNEL);
if strcmp(fPA,'firstPhaseAmplitude=-1') || strcmp(fPD,'firstPhaseDuration=-1') % stop the trial
    fprintf(t,'enabled=0');
else
    fprintf(t,'enabled=1');
end
fprintf(t,pTD);
fprintf(t,pT);
fprintf(t,trEL);
fprintf(t,nPT);
fprintf(t,pTP);
fprintf(t,pSR);
fprintf(t,fPD);
fprintf(t,sPD);
fprintf(t,IPD);
fprintf(t,fPA);
fprintf(t,sPA);
% fprintf(t,stS);
fprintf(t,'stimShape=1'); % Stimulation Shape = 1, Biphasic + Inter-phase delay
fprintf(t,stP);
fprintf(t,enAS);
fprintf(t,prAS);
fprintf(t,poAS);
fprintf(t,maAS);
fprintf(t,enCR);
fprintf(t,prCR);
fprintf(t,poCR);
fprintf(t,['ID=' num2str(n)]);
end


