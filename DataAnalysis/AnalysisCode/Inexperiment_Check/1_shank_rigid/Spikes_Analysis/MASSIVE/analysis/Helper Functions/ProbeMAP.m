%% Generates a NeuroProbe MAP variable
function E_MAP = ProbeMAP
% Initialise
E_MAP = cell(33,1);
E_MAP(1,1) = {'ATLAS'};
E_MAP(1,2) = {'CBC-Series NN'};
E_MAP(1,3) = {'I-Series NN'};
E_MAP(1,4)={'2 Shank I-Series NN'};
E_MAP(1,5)={'FLEXIBLE I-Series NN'};

MOLC_MALE = [32,30,31,28,29,27,25,22,23,21,17,18,19,20,24,26,1,4,13,14,15,16,12,10,8,6,2,3,5,7,9,11];
MOLC_FEMALE = [32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1];
OMNETICS_MALE = [23,25,27,29,31,19,17,21,11,15,13,1,3,5,7,9,10,8,6,4,2,14,16,12,22,18,20,32,30,28,26,24];
INTAN = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
ATLAS = [28,12,30,14,26,10,03,19,05,21,31,15,29,13,01,17,00,16,27,11,02,18,07,23,04,20,25,09,24,08,06,22];
NEURONEXUS_A = [22,9,29,3,14,2,12,6,17,30,20,25,15,31,1,16,0,21,5,10,26,23,7,8,24,11,28,13,4,19,18,27];
NEURONEXUS_I = [22,27,18,28,13,4,9,29,19,3,14,2,12,6,17,30,20,25,15,31,11,1,16,0,21,5,10,26,23,7,8,24];
NN_I = [1,32,2,31,3,30,4,29,5,28,6,27,7,26,8,25,9,24,10,23,11,22,12,21,13,20,14,19,15,18,16,17];
NN_I2shank=[22,8,18,23,13,10,9,21,19,16,14,11,12,15,17,20,24,27,7,28,26,4,5,29,0,3,1,2,31,6,25,30];
NN_flexible=[8,7,9,6,10,5,12,3,13,2,14,1,23,24,22,25,21,26,19,28,18,29,17,30,16,31,20,27,15,0,11,4];
for n = 1:32
    if (ATLAS(n) < 10)
        temp = strcat('A-00',num2str(ATLAS(n)));
    else
        temp = strcat('A-0',num2str(ATLAS(n)));
    end
    E_MAP(n+1,1) = {temp};   
    if (NEURONEXUS_A(n) < 10)
        temp = strcat('A-00',num2str(NEURONEXUS_A(n)));
    else
        temp = strcat('A-0',num2str(NEURONEXUS_A(n)));
    end
    E_MAP(n+1,2) = {temp};
    M = find(INTAN == find(OMNETICS_MALE == find(MOLC_FEMALE == find(MOLC_MALE == NN_I(n)))));
    if (M-1 < 10)
        temp = strcat('A-00',num2str(M-1));
    else
        temp = strcat('A-0',num2str(M-1));
    end
    E_MAP(n+1,3) = {temp};
    L=NN_I2shank(n);
    if (L < 10)
        temp = strcat('A-00',num2str(L));
    else
        temp = strcat('A-0',num2str(L));
    end
    E_MAP(n+1,4) = {temp};
    T=NN_flexible(n);
    if (T < 10)
        temp = strcat('A-00',num2str(T));
    else
        temp = strcat('A-0',num2str(T));
    end
    E_MAP(n+1,5) = {temp};
    
end
end