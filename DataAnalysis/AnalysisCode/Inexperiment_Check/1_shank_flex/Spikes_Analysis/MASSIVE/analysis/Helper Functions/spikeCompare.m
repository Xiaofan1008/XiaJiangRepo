function chk = spikeCompare(wave,box)
%% Checks to see if a spike waveform passes through a series of boxes
Xs = 1:1:49;
chk = zeros(size(wave,1),1);
%% Logic
% Alignment
[~,m] = min(wave,[],2);
chk(m == 13) = 1;
% Boxing
tmpbox1 = box(1,1:5);
tmpbox2 = box(1,6:10);
for w = 1:size(wave,1)
    try [xi,~] = polyxpoly(Xs,wave(w,:),tmpbox1,tmpbox2);
    catch
        continue;
    end
    for i = 1:size(box,1)        
        [xi,~] = polyxpoly(Xs,wave(w,:),box(i,1:5),box(i,6:10));
        if isempty(xi)
            chk(w) = 0;
            continue;
        end
    end
end
end