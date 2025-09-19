function data = complexBlank(data)
%% Variables
threshfac = 1.2;
warning off signal:findpeaks:largeMinPeakHeight
%% This function is designed to locate abrupt changes in potential and smooth the waveform out to prevent filter ringing.
for c = 1:size(data,1)
    thresh = std(data(c,:))*threshfac;
    pt = abs(diff(data(c,:)));
    while ~isempty(pt)
        pt = abs(diff(data(c,:)));
        [~,pt] = findpeaks(pt,'MinPeakHeight',thresh,'MinPeakDistance',3000);
        for i = 1:length(pt)
            thresh1 = interpolate(data(c,pt(i):pt(i)+3000),1);
            thresh2 = interpolate(data(c,pt(i)+1:pt(i)+3001),1);
            data(c,pt(i):pt(i)+3000) = data(c,pt(i)+1:pt(i)+3001) - thresh2 + thresh1;
        end
    end
end