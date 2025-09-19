% This script is used to verify spike detection
dbstop on error
clear all
%% Variables and Loading
loadTrig;
loadSpikes;
trig = cast(trig,'int32');
n = 20;
BN = [-20 20];
nChn = 32;
FS = 30000;
nT = length(trig);
T = randi(nT,1,n);
XM = 1/30:1/30:(diff(BN)*n);
d = Depth(3);
%% Generate the mu vector
mu_dir = dir('*.mu.dat');
mu_name = mu_dir.name;
mu_bytes = mu_dir.bytes;
L = mu_bytes / (2*FS*nChn);
t = 60;
mu_fid = fopen(mu_name);
figure; hold on;
for c = 2:5
    disp(['Running channel ' num2str(c)]);
    mu = zeros(1,L*FS);
    fseek(mu_fid,0,'bof');
    for l = 1:ceil(L/t)
        v = fread(mu_fid,[nChn t*FS],'short') ./ 10;
        try
            mu(1,1+((l-1)*t*FS):((l)*t*FS)) = v(d(c),:);
        catch
            mu(1,1+((l-1)*t*FS):end) = v(d(c),:);
        end
    end
    %% Logic
    SPIKES = cell(1,n);
    spk = sp{d(c)};
    %spk = denoiseSpikes(spk);
    spk = spk(:,1);
    for i = 1:n
        M(1,1+(i-1)*(diff(BN)*(FS/1e3)):(i)*(diff(BN)*(FS/1e3))) = mu(1,1+(trig(T(i))+BN(1))*(FS/1e3):(trig(T(i))+BN(2))*(FS/1e3));
        SP = spk(spk > trig(T(i))+BN(1) & spk < trig(T(i))+BN(2));
        if ~isempty(SP)
            for k = 1:length(SP)
                SP(k) = SP(k) - trig(T(i));
            end
        end
        SPIKES{i} = SP;
    end
    %% Plotting
    plot(XM,M + 200*(c-1),'Color','k');
    for i = 1:n
        spk = cell2mat(SPIKES(i));
        line([abs(BN(1))+(i-1)*diff(BN) abs(BN(1))+(i-1)*diff(BN)],[-200 + 200*(c-1) 200 + 200*(c-1)],'Color','b');
        for s = 1:length(spk)
            text(spk(s)+(abs(BN(1))+diff(BN)*(i-1)),60+200*(c-1),'*','FontSize',16);
        end
    end
end
%% Close