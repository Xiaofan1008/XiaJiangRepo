%% This function takes a waveform and calculates instantaneous phase along it
function phase = calcPhase(wave,band)
% User-defined variables
SPECT = 0;
PSD = 0;
FS = 1000;
tapers = [0.05,20];
if (band == 4)
    tapers = [0.5,2];
end
if (band == 5)
    tapers = [0.2,5];
end
if (band == 6)
    tapers = [0.5,2];
end
if (band == 10)
    tapers = [0.2,5];
end
if (band == 20)
    tapers = [0.1,10];
end
if (band == 80)
    tapers = [0.025,40];
end
%% Butterworth filter for debugging
%[b,a] = butter(10,[40 120]/500);
%% Taper explanation:
% tapers(2) is the bandwidth of the filter. The filter will allow
% frequencies from 'band' +/- tapers(2).
% tapers(1) is the length of data we're giving the filter. The filter
% window is equal to tapers(1) * FS;
% What is critical when using this function is that the stimulation event
% is not included in the waveform passed in - wave should stop 1 sample
% before the stimulation event.
%% Use mtfilter to generate a bandpass filtered waveform with complex values
phase = wave;
for i = 1:size(wave,1)
    Y = mtfilter(wave(i,:),tapers,FS,band,0,1);
    %Z = fliplr(filter(b,a,fliplr(wave(i,:))));
    phase(i,:) = angle(Y).*(180/pi); % Angle returns radians, by default
end
if (PSD)
    disp('Hi');
    Fs = 1000;
    t = 0:1/Fs:1-1/Fs;
    xdft = fft(Y);
    N = length(t);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(t):Fs/2;
    plot(freq,10*log10(psdx))
end
% We also want to be able to output a spectrogram of the LFP at each
% trigger line. This allows us to see how isolated each frequency band is
% from artefact
if (SPECT)
    window = [-500 500]; trig = []; loadTrig;
    TrialParams = []; trialID = 1; loadTrialParams;
    TrialParams = TrialParams(2:end,:);
    TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == trialID));
    trig = trig(TrialParams);
    trig = cast(trig,'int32');
    nTrig = length(trig);
    x = zeros(1,diff(window)+1);
    for t = 1:nTrig
        x = x + real(Y(trig(t)+window(1):trig(t)+window(2)));
    end
    x = x ./ t;
    spectrogram(x,32,4,128,1e3,'yaxis');
end
end