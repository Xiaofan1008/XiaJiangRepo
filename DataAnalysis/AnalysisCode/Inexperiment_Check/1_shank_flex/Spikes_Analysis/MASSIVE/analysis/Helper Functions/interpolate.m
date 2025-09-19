function wave = interpolate(wave,order)
XSTART = 1; XSTOP = length(wave);
X = [XSTART,XSTOP];
[CF, S, mu] = polyfit(X,[wave(XSTART),wave(XSTOP)],order);
X = XSTART:1:XSTOP;
Y = polyval(CF,X,S,mu);
wave(XSTART:XSTOP) = Y;
end