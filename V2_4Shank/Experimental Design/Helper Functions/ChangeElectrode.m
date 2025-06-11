function ChangeElectrode()

home = pwd; %returns the current directory in the string "home"
[DATAFILE,FOLDER] = uigetfile; % The experimental protocol file
cd(FOLDER);
oldName = input('Please enter the old channel name like so: PAS1E5,PBS2E10,...,PDS3E13 \nWhere P is port letter, S is #shank, E is #electrode\n','s'); % channels on port B are +32 shanks are 16 electrodes each
newName = input('Please enter the new channel name like so: PAS1E5,PBS2E10,...,PDS3E13 \nWhere P is port letter, S is #shank, E is #electrode\n','s'); % channels on port B are +32 shanks are 16 electrodes each
% Local variables
exp_datafile = dir(DATAFILE);

% Load in the stimulation parameters. 
data = load(exp_datafile.name);
StimParams = data.StimParams;

% Check header to ensure channel name is in the first column
header = StimParams(1,:);
if ~strcmp(header{1}, "CHANNEL")
    error('Expected CHANNEL to be in the first column');
end

% Initialize
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,data.E_Mapnumber+5);

% MAP oldName to electrode
port_shank_electrode=split(oldName,["P","S","E",","]);
shank_electrode=(cellfun(@(x) str2double(x),port_shank_electrode,'UniformOutput',false));
Portnum=(cellfun(@(x) x(isletter(x)),port_shank_electrode,'UniformOutput',false));
Portnum=Portnum(~cellfun('isempty',Portnum));
shank_electrode=cell2mat(shank_electrode);
shank_electrode=shank_electrode(~isnan(shank_electrode));
loopcounter=0;
CHN=zeros(1,length(shank_electrode)/2);
for shank=1:2:length(shank_electrode)
    loopcounter=loopcounter+1;
    CHN(loopcounter)=shank_electrode(shank+1);
    if shank_electrode(shank)==4 || shank_electrode(shank)==3
        CHN(loopcounter)=CHN(loopcounter)+16;
    elseif shank_electrode(shank)>4
        error('Shank does not exist')
    end
    if strcmp(Portnum(loopcounter),'B')
        CHN(loopcounter)=CHN(loopcounter)+32;
    elseif strcmp(Portnum(loopcounter),'C')
        CHN(loopcounter)=CHN(loopcounter)+64;
    elseif strcmp(Portnum(loopcounter),'D')
        CHN(loopcounter)=CHN(loopcounter)+96;
    end
end

oldName = E_MAP{CHN+1,1};

% MAP newName to electrode
port_shank_electrode=split(newName,["P","S","E",","]);
shank_electrode=(cellfun(@(x) str2double(x),port_shank_electrode,'UniformOutput',false));
Portnum=(cellfun(@(x) x(isletter(x)),port_shank_electrode,'UniformOutput',false));
Portnum=Portnum(~cellfun('isempty',Portnum));
shank_electrode=cell2mat(shank_electrode);
shank_electrode=shank_electrode(~isnan(shank_electrode));
loopcounter=0;
CHN=zeros(1,length(shank_electrode)/2);
for shank=1:2:length(shank_electrode)
    loopcounter=loopcounter+1;
    CHN(loopcounter)=shank_electrode(shank+1);
    if shank_electrode(shank)==4 || shank_electrode(shank)==3
        CHN(loopcounter)=CHN(loopcounter)+16;
    elseif shank_electrode(shank)>4
        error('Shank does not exist')
    end
    if strcmp(Portnum(loopcounter),'B')
        CHN(loopcounter)=CHN(loopcounter)+32;
    elseif strcmp(Portnum(loopcounter),'C')
        CHN(loopcounter)=CHN(loopcounter)+64;
    elseif strcmp(Portnum(loopcounter),'D')
        CHN(loopcounter)=CHN(loopcounter)+96;
    end
end

newName = E_MAP{CHN+1,1};

% Replace oldName with newName
for i = 2:size(StimParams, 1)
    if strcmp(StimParams{i, 1}, oldName)
        StimParams{i, 1} = {newName};
    end
end

data.StimParams = StimParams;

% Save
here = pwd;
fileID = '001';
DIR = uigetdir;
C = strsplit(DIR,'\');
NAME = C{size(C,2)};
NEWFILE = strcat(NAME,'_exp_datafile_',fileID,'.mat');
id = 1;
cd(DIR);
dataPlace = pwd;
save(NEWFILE, '-struct', 'data');
disp(['Experimental datafile: ' NAME '_exp_datafile_Changed_Electrode_' fileID '.mat has been saved']);
cd (here);

end

