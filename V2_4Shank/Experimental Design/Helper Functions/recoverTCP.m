function recoverTCP(t)
%% A recovery mode function that clear and stabilises the connection
CHK = 0;
while(~CHK)
    flushinput(t)
    fprintf(t,'RECOVER');
    pause(1);
    data = fscanf(t);
    data = string(data(1:end-1));
    if strcmp(data,"PHEW")
        CHK = 1;    
    end
    flushinput(t)
end
end