function y = ste(x)
% Compute standard error of given sample distribution x
STD = std(x);
n = length(x);
y = STD/sqrt(n);
end