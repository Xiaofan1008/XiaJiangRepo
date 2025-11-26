function depth = Depth(varargin)
filepath = pwd;
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
NN_cutoff = datetime('06-Nov-2018 00:00:00');
Flex_cutoff = datetime('14-Jul-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);

if datetime(fileinfo.date) < NN_cutoff
    NN_order = 0;
else
    NN_order = 1;
end
if datetime(fileinfo.date) > datetime('28-Nov-2018 00:00:00')
    NN_I = 1;
else 
    NN_I = 0;
end
if ~(NN_order)
    p = 1;
else
    p = 2;
end
if (NN_I)
    p = 3;
end
if nChn>=64
    p=6;
elseif datetime(fileinfo.date) > Flex_cutoff
    p = 5;
end
if nargin==1
    E_Mapnumber=cell2mat(varargin(1));
    if E_Mapnumber>0 || nChn>64
        p=6;
    elseif datetime(fileinfo.date) > Flex_cutoff
        p=5;
    else
        p=3;
    end
end
E_MAP = ProbeMAP;
depth = zeros(1,nChn);
for n = 2:nChn+1
    str = E_MAP{n,p};
    if strcmp(str(1),'A')
        str = str(4:5);
        str = str2double(str);
    elseif strcmp(str(1),'B')
        str = str(4:5);
        str = str2double(str)+32;
    elseif strcmp(str(1),'C')
        str = str(4:5);
        str = str2double(str)+64;
    elseif strcmp(str(1),'D')
        str = str(4:5);
        str = str2double(str)+96;
    end
    depth(n-1) = str;
end
depth = depth';
if ~isempty(depth == 0)
    depth = depth + 1;
end
end