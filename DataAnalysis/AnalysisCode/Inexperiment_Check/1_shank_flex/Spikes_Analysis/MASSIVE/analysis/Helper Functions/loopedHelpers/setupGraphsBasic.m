if ~exist('nChn','var')
    nChn = 32;
end
if nChn(end) == 32
    ROW = 6;
    COL = 6;
    A = zeros(6,6);
    for i = 1:32
        A(i+3) = i;
    end
    for i = 1:28
        A(i+2) = i;
    end
    A(1,6) = 0;
    for i = 1:4
        A(i+1) = i;
    end
    A(6,1) = 0;
end
XLABEL = [25,30,32,33,34,35];
YLABEL = [7,19,32];