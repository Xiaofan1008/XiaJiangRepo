%% Script that repairs broken experimental datafiles
path = uigetdir;
RECORDINGS = dir(path);
RECORDINGS(1:2) = [];
nDATA = zeros(1,length(RECORDINGS));
TrialParams = [];
E_MAP = AtlasMAP;
C = E_MAP(2:end,5);
for i = 1:length(RECORDINGS)
    datafile = dir([path '\' RECORDINGS(i).name '\*_exp_datafile*.mat']);
    if ~isempty(datafile)
        load([path '\' RECORDINGS(i).name '\' datafile.name]);
        if exist('TrialParams','var')
            clear TrialParams;
            continue;
        else
            TrialParams = cell(length(StimParams(:,1)),3);
            TrialParams{1,1} = 'Trial Number';
            TrialParams{1,2} = 'Trial ID';
            TrialParams{1,3} = 'Channel';
            for n = 2:length(StimParams(:,1))
                TrialParams{n,1} = n-1;
                TrialParams{n,2} = (find(AMP == StimParams{n,16},1)) + (length(AMP)*(find(DUR == StimParams{n,13},1)-1)) + ...
                (length(AMP)*length(DUR)*(find(CHN == find(ismember(C,StimParams{n,1}),1),1)-1));
                TrialParams{n,3} = StimParams{n,1};
            end
            save([path '\' RECORDINGS(i).name '\' datafile.name],'TrialParams','StimParams','n_Trials','E_MAP','E_MAT','E_DIAM','n_REP','rand_order','AMP','DUR','CHN','PARAM','PULSE');
        end
    end
end