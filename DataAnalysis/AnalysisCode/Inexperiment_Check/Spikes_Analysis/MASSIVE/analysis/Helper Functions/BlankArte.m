%% Function analyses a provided waveform and blanks the artifact
function blankedwave = BlankArte(wave,BIN,count)
[XSTART,XSTOP] = BlankingHelper(wave,BIN);
if isempty(XSTART) && isempty(XSTOP)
    blankedwave = wave;
    return;
end
XSTART = floor(XSTART*30 - 0.5);
XSTOP = ceil(XSTOP*30 + 1);
X = [XSTART,XSTOP];
[CF, S, mu] = polyfit(X,[wave(XSTART),wave(XSTOP)],1);
X = XSTART:1:XSTOP;
Y = polyval(CF,X,S,mu);
blankedwave = wave;
try
    blankedwave(XSTART:XSTOP) = Y;
catch
    blankedwave(XSTART:XSTOP) = Y';
end
if max(abs(blankedwave)) > 1e3 % Try again.
    if (count > 5)
        % Give up
        return;
    end
    blankedwave = BlankArte(blankedwave,BIN,count+1);
end
end