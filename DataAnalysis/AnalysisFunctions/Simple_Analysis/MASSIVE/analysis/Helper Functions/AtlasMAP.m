%% Generates an Atlas NeuroProbe 208 MAP file
function E_MAP = AtlasMAP
% Initialise
E_MAP = cell(33,4);
E_MAP(1,1) = {'Electrode'};
E_MAP(1,2) = {'MOLC_MALE'};
E_MAP(1,3) = {'MOLC_FEMALE'};
E_MAP(1,4) = {'OMNETICS MALE'};
E_MAP(1,5) = {'INTAN'};

ELECTRODES = (1:1:32);
MOLC_MALE = [3,33,7,37,18,28,4,34,17,27,9,39,5,35,8,38,10,40,1,31,6,36,19,29,2,32,16,26,20,30,15,25];
MOLC_FEMALE = [30,3,26,7,20,13,29,4,19,14,24,9,28,5,25,8,23,10,32,1,27,6,21,12,31,2,18,15,22,11,17,16];
OMNETICS_MALE = [19,3,17,1,21,5,12,28,10,26,16,0,18,2,14,30,15,31,20,4,13,29,8,24,11,27,22,6,23,7,9,25];
INTAN = [28,12,30,14,26,10,03,19,05,21,31,15,29,13,01,17,00,16,27,11,02,18,07,23,04,20,25,09,24,08,06,22];
for n = 1:32
    E_MAP(n+1,1) = {ELECTRODES(n)};
    E_MAP(n+1,2) = {MOLC_MALE(n)};
    E_MAP(n+1,3) = {MOLC_FEMALE(n)};
    E_MAP(n+1,4) = {OMNETICS_MALE(n)};
    if (INTAN(n) < 10)
        temp = strcat('A-00',num2str(INTAN(n)));
    else
        temp = strcat('A-0',num2str(INTAN(n)));
    end
    E_MAP(n+1,5) = {temp};
end
end