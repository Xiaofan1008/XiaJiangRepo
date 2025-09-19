function t = subsample(t,n)
%% Function allows subsampling randomly from a vector
rng('shuffle');
l = length(t);
tl = randperm(l,n);
t = t(tl);
end