function lfp = lfpBlank(lfp,trig)
%% This function takes an lfp waveform and set of trigger lines and blanks the artefact resulting from stimulation triggers
%trig = cast(trig,'int32');
%% Blanking variables
blankWINDOW = [-1 21];
%% Blanking procedures
for t = 1:length(trig)
    XSTART = trig(t) + blankWINDOW(1);
    XSTOP = trig(t) + blankWINDOW(2);
    YSTART = lfp(cast(XSTART,'int32'));
    YSTOP = lfp(cast(XSTOP,'int32'));
    [CF, S, mu] = polyfit([XSTART,XSTOP],[YSTART,YSTOP],1);
    X = XSTART:1:XSTOP;
    Y = polyval(CF,X,S,mu);
    lfp(cast(XSTART,'int32'):cast(XSTOP,'int32')) = Y;
end
end