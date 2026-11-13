function OK = checkTCPOK(t)
% Initialise the check
flushinput(t)
fprintf(t,'RUOK');
data = fscanf(t);
data = string(data(1:end-1));
if ~strcmp(data,"YEAH")
    OK = 0;
else
    OK = 1;
end
flushinput(t)
end
