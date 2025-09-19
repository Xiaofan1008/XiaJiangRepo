% Helper script for assessing waveforms derived from phy
dbstop on error
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*mu.dat']);
if ~isempty(tmp)    
    sz = tmp.bytes/2;     
    tmp = [filepath '\' tmp.name];
    m_fid = fopen(tmp,'r');
    nLoops = ceil(sz / ((30*1440000)/32));
    wave = zeros(size(sp,1),51);
    for n = 1:nLoops
        mu = fread(m_fid,[32, (30*1440000)/32],'short');
        if isempty(mu)
            return;
        end
        mu = mu ./ 10;
        spikes = cast(sp(sp < size(mu,2)*(n) & sp > size(mu,2)*(n - 1)),'int32');
        spikes = spikes - size(mu,2)*(n-1);
        if spikes(end) > (1350000-25)
            spikes(end) = [];
        end
        if spikes(1) < 26
            spikes(1) = [];
        end
        nWaves = length(spikes);
        M = mean(mu(29,:));
        for i = 1:nWaves
            wave(i,:) = mu(29,spikes(i)-25:spikes(i)+25) - M;
        end
    end
end
