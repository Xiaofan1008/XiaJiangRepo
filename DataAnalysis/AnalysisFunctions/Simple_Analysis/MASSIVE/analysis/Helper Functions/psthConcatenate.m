function psthConcatenate
dbstop if error
%% Variables
BIN = [-15 50]; % ms
chan = 6;
SPACING = 0;
tID = 8;
phaseBins = [360 6];
%% Loading in data
d = Depth; sp = loadSpikes; sp = sp{d(chan)};
figure; hold on;
for tt = tID
    trig = loadTrig(0);
    output = trig_helper([6 -251],[6]);
    psthTrig = cell(1,phaseBins(2));
    for p = 1:phaseBins(2)
        ind = output((output(:,5) == p),:);
        ind = ind(ind(:,2) == tt,:);
        [~,chk] = unique(round(ind(:,3)));
        ind = ind(chk,:);
        psthTrig{p} = trig(ind(:,1));
        nTrig = length(psthTrig{p});
        clear SPIKES; SPIKES = cell(1,nTrig);
        for t = 1:nTrig
            SPIKES{t} = sp(sp(:,1) >= psthTrig{p}(t)/30 + BIN(1) & sp(:,1) <= psthTrig{p}(t)/30 + BIN(2),1) - psthTrig{p}(t)/30 - BIN(1);
        end
        if p==1; psth(SPIKES,BIN,2,nTrig,[],0,[0 SPACING]); SPACING = SPACING + nTrig;
        else; psth(SPIKES,BIN,2,nTrig,[],0,[0 SPACING+4]); SPACING = SPACING + nTrig + 4; end
        line([BIN(1) BIN(2)],[SPACING+2 SPACING+2],'Color','b');
    end
    line([0 0],[0 SPACING+2],'Color','r');
end
end