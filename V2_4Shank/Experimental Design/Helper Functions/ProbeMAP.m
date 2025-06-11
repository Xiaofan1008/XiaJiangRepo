%% Generates a NeuroProbe MAP variable
function E_MAP = ProbeMAP
% Initialise
E_MAP = cell(33,1);
E_MAP(1,1) = {'ATLAS'};
E_MAP(1,2) = {'CBC-Series NN'};
E_MAP(1,3) = {'I-Series NN'};
E_MAP(1,4)={'2 Shank I-Series NN'};
E_MAP(1,5)={'FLEXIBLE I-Series NN'};
E_MAP(1,6)={'4SHANK FLEXIBLE I-Series NN'};


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

NN_flexible_4SHANKA=[31,7,0,24,30,6,1,25,29,5,2,26,28,4,3,27,  8,16,23,15,9,17,22,14,10,18,21,13,11,19,20,12]; %S1 and S4
NN_flexible_4SHANKB=[8,11,9,15,10,19,12,23,13,22,14,21,16,20,17,18,   4,7,0,6,28,5,24,3,25,2,26,1,27,31,29,30]; %S2 and S3
NN_flexible_4SHANKC=[31,7,0,24,30,6,1,25,29,5,2,26,28,4,3,27,  8,16,23,15,9,17,22,14,10,18,21,13,11,19,20,12]; %S1 and S4
NN_flexible_4SHANKD=[8,11,9,15,10,19,12,23,13,22,14,21,16,20,17,18,   4,7,0,6,28,5,24,3,25,2,26,1,27,31,29,30]; %S2 and S3


%ECOG
omneticsconnector=1:32;
nnconnector=fliplr([1 3 5 7 9 11 13 15 50 52 54 56 58 60 62 64 63 61 59 57 55 53 51 49 16 14 12 10 8 6 4 2]);
chnsordershank=[1 16 2 15 3 14 4 13 5 12 6 11 7 10 8 9, 49 64 50 63 51 62 52 61 53 60 54 59 55 58 56 57];
nnconnector=[34 43 44 45 33 46 47 48 17 18 19 32 20 21 22 31 23 24 25 30 26 27 28 29 36 37 38 39 35 40 41 42];
chnsordershank=[17 32 18 31 19 30 20 29 21 28 22 27 23 26 24 25, 33 48 34 47 35 46 36 45 37 44 38 43 39 42 40 41];
order=zeros(32,1);
for i=1:32
    order(i)=find(chnsordershank(i)==nnconnector);
end


ECOG_order=[1 15 2 14 3 13 4 12 5 11 6 10 7 9 8,   16 29 17 28 18 27 19 26 20 25 21 24 22 23,   30 43 31 42 32 41 33 40 34 39 35 38 36 37,   44 57 45 56 46 55 47 54 48 53 49 52 50 51,   58 71 59 70 60 69 61 68 62 67 63 66 64 65,  ...
    72 85 73 84 74 83 75 82 76 81 77 80 78 79,    86 99 87 98 88 97 89 96 90 95 91 94 92 93,    100 113 101 112 102 111 103 110 104 109 105 108 106 107,    114 128 115 127 116 126 117 125 118 124 119 123 120 122 121];
nnconnector=[63 50 47 34 31 18 15 2 49 64 33 48 17 32 1 16 14 3 30 19 46 35 62 51 4 13 20 29 36 45 52 61,    8 9 24 25 40 41 56 57 10 7 26 23 42 39 58 55 53 60 37 44 21 28 5 12 59 54 43 38 27 22 11 6,...
     128 113 112 97 96 81 80 65 114 127 98 111 82 95 66 79 77 68 93 84 109 100 125 116 67 78 83 94 99 110 115 126,   69 76 85 92 101 108 117 124 75 70 91 86 107 102 123 118 120 121 104 105 88 89 72 73 122 119 106 103 90 87 74 71];
order=zeros(128,1);
 for i=1:128
    order(i)=find(ECOG_order(i)==nnconnector);
end

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
for n = 1:128
    if n<33
        T=NN_flexible_4SHANKA(n);
        if (T < 10)
            temp = strcat('A-00',num2str(T));
        else
            temp = strcat('A-0',num2str(T));
        end
    elseif n<65
        T=NN_flexible_4SHANKB(n-32);
        if (T < (10))
            temp = strcat('B-00',num2str(T));
        else
            temp = strcat('B-0',num2str(T));
        end
    elseif n<97
        T=NN_flexible_4SHANKC(n-64);
        if (T < (10))
            temp = strcat('C-00',num2str(T));
        else
            temp = strcat('C-0',num2str(T));
        end
    else
        T=NN_flexible_4SHANKD(n-96);
        if (T < (10))
            temp = strcat('D-00',num2str(T));
        else
            temp = strcat('D-0',num2str(T));
        end
    end
    E_MAP(n+1,6) = {temp};
    
end 
end