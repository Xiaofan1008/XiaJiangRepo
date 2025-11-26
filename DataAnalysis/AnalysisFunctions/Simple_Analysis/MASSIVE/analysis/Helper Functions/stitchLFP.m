function stitchLFP
%% This function generates an amplifier_lfp.dat file where stimulation and parameter updates
% are stitched over
file = pwd;
%% Variables
FS = 30000;
nChn = 32;
T = 65;
BIN = [-10, 250] .* 30;
OVERLAP = 3e3;
%% Load in data
v_fid = fopen([file '\amplifier.dat'],'r');
lfp_fid = fopen([file '\amplifier_lfp.dat'],'w');
loadTrig;
v_dir = dir([file '\amplifier.dat']);
bytes = v_dir.bytes;
LEN = bytes / 2 / nChn / FS;
nL = ceil(LEN / T);
%% Grab a bit of data to use
v = fread(v_fid,[nChn, 5*FS],'int16') .* 0.195;
stitch = zeros(nChn,diff(BIN)+1);
for c = 1:nChn
    for w = 1:50:diff(BIN)-51
        if (range(v(c,w:w+diff(BIN))) < 400)
            stitch(c,:) = v(c,w:w+diff(BIN));
            break;
        end
    end
end
chk(:) = sum(stitch(:,:),2);
for c = 1:nChn
    if sum(stitch(c,:)) == 0
        stitch(c,:) = stitch(find(chk ~= 0,1),:);
    end
end
%% Stitch data
fseek(v_fid,0,'bof');
for n = 1:nL
    disp(['Running over loop ' num2str(n)]);
    thisTrig = trig(trig > ((n-1)*T*1e3)+(OVERLAP/3) & trig < n*T*1e3-(OVERLAP/3)) - ((n-1)*T*1e3);
    thisTrig = cast(thisTrig,'int32');
    nTrig = length(thisTrig);
    offset = ((n-1)*T*FS*nChn*2);
    fseek(v_fid,offset,'bof');
    v = fread(v_fid,[nChn, (T+OVERLAP/1e3)*FS],'int16') .* 0.195;
    if isempty(v)
        break;
    end
    for c = 1:nChn
        for t = 1:nTrig            
            v(c,(thisTrig(t)*30+BIN(1)):(thisTrig(t)*30+BIN(2))) = stitch(c,:) + diff([stitch(c,1),v(c,(thisTrig(t)+BIN(1))*30 - 1)]);
        end
    end
    fwrite(lfp_fid,v ./ 0.195,'int16');            
end
%% Close down
fclose(v_fid);
fclose(lfp_fid);
end