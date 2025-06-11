%% Calculates an approximation of the actual waveform delivered by Intan, given stimulus parameters
function error = Intan_error(amp,dur,R,gap,tau)
% Calculate the desired square wave
x = linspace(1,1000,1000);
y = zeros(1,1000);

% Calculate the artifact
y(100) = -9;
y(100+dur+gap) = 9;

for n = 101:100+dur
    y(n) = -amp*R + (-9 - -amp*R)*(exp(-(n-100)/tau));
end
for n = 101+dur:101+dur+gap
    y(n) = 0 + (y(100+dur) - 0) * exp(-(n-(101+dur))/tau);
end
for n = 101+dur+gap:101+dur+gap+dur
    y(n) = amp*R + (9 - amp*R) * exp(-(n-(101+dur+gap))/tau);
end
for n = 101+gap+(2*dur)+1:101+gap+dur+gap+dur
    y(n) = 0 + (y(100+(2*dur)+gap) - 0) * exp(-(n-(101+(2*dur)+gap))/tau);
end
plot(x,y)
z = trapz(x,y);
error = y;
end