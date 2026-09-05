function v = artefactBlank(v)
dbstop on error
%% User Variables
ARTMAX = 5e2;
%% SETUP
for c = 1:size(v,1)
    % Remove DC offset
    v(c,:) = v(c,:) - mean(v(c,:));
    % Count the number of artefacts
    indexArt = find(abs(v(c,:)) >= ARTMAX);
    chk = find(diff(indexArt) > 750);
    %chk = [chk length(indexArt)];
    tmp = zeros(1,length(chk));
    for i = 1:length(chk)
        iStart = 1;
        iStop = chk(i);
        if (i > 1)
            iStart = chk(i-1)+1;
        end        
        tmp(i) = find(abs(v(c,indexArt(iStart):indexArt(iStop))) == max(abs(v(c,indexArt(iStart):indexArt(iStop)))),1);
        tmp(i) = ceil((tmp(i) / (indexArt(iStop) - indexArt(iStart))) * (iStop - iStart) + iStart);
    end
    tmp(isnan(tmp)) = [];
    indexArt = indexArt(tmp);
    nArt = length(indexArt) + 1;
    while (nArt ~= 1)
        nArt = nArt - 1;
        % For each artefact, calculate an appropiate window
        WINDOW = [0 0];
        if indexArt(nArt) > 1500
            WINDOW(1) = indexArt(nArt) - 1500;
        else
            continue;
        end
        if indexArt(nArt) < (length(v(c,:)) - 1500)
            WINDOW(2) = indexArt(nArt) + 1500;
        else
            continue;
        end
        BIN = [-(diff(WINDOW)/60) (diff(WINDOW)/60)];
        % For each artefact, call blank        
        v(c,WINDOW(1):WINDOW(2)) = BlankArte(v(c,WINDOW(1):WINDOW(2)),BIN,1);     
    end    
end

end

