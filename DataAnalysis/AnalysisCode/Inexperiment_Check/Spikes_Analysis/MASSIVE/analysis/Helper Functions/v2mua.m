[Mufilt,~] = generate_Filters;
MuNf = length(Mufilt);
nData = size(v,2);
mu = v;
for c = 1:size(v,1)
    tmp = conv(v(c,:),Mufilt);
    mu(c,:) = tmp(1,MuNf/2:nData+MuNf/2-1);
end
clear Mufilt MuNf nData tmp c