if ~isempty(dir('*exp_datafile*.mat'))
    Tname = dir('*exp_datafile*.mat');
    Tname = Tname.name;
    StimParams = load(Tname,'StimParams');
    StimParams = StimParams.StimParams;
    delay = StimParams(2,19);
    delay = delay{1};
else
    delay = 0;
end