function S = verifyParameters(t,n)
S = 0;
data = fscanf(t);
data = string(data(1:end-1));
if ~strcmp(data,strcat("Parameter ", num2str(n), " Set"))
    disp('ERROR. Settings not yet set.')
else
    S = 1;
end
flushinput(t);
end