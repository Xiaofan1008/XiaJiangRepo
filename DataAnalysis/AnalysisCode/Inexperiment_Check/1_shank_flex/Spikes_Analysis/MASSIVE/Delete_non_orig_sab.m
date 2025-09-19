function Delete_non_orig_sab
%Deletes all non-original files from file path

fclose('all');
if ~isempty(dir('*_dn_sab.dat'))
    delete *_dn_sab.dat
end

if ~isempty(dir('*.mu_sab.dat'))
    delete *.mu_sab.dat
end

if ~isempty(dir('*.sp.mat'))
    delete *.sp.mat
end

if ~isempty(dir('*trig.dat'))
    delete *trig.dat
end
end

