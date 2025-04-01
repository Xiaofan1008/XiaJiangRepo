% Helper script for loading trigger times
function trig = loadTrig(chk)

    global DATA_FOLDER;
    tmp = dir(fullfile(DATA_FOLDER, '*.trig.dat'));
    if isempty(tmp)
       error('No trigger file found in %s', DATA_FOLDER);
    end

    trig_file = fullfile(DATA_FOLDER, tmp(1).name);
    t_fid = fopen(trig_file, 'r');
    trig = fread(t_fid, [1, inf], 'double');
    fclose(t_fid);
end



% if ~nargin
%     chk = 0;
% end
% filepath = pwd;
% tmp = dir([filepath filesep '*.trig.dat']);
% if ~isempty(tmp)
%     sz = tmp.bytes/2;
%     tmp = [filepath filesep tmp.name];
%     t_fid = fopen(tmp,'r');
%  trig = fread(t_fid,[1, sz],'double'); % byte rate = 18560; 'double' = 64bits = 8 bytes; 
%     fclose(t_fid);
% else    
%     return;
% end
% 
% if isempty(trig)
%     return
% end
% if ~(chk)
%     return;
% elseif (chk==1)
%     trig = trig ./ 30;
%     trig = cast(trig,'int64');
% end
% end