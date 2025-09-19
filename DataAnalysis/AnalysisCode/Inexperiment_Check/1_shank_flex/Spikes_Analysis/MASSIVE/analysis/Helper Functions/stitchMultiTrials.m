function stitchMultiTrials
%% This function is passed multiple Intan and/or Blackrock datasets and experimental datafiles
% and then stitches them together such that all trials are present in a
% single directory
basePath = 'X:\Tim\';
%% Datasets
% Intan
IntanlOD = {'CJ202\estim_pen7_003_191004_133556\',...
            'CJ202\estim_pen7_003_191004_144648\',...
           };
% Blackrock
BlackrocklOD = {'CJ202\estim_pen7_003_191004_133556\',...
                'CJ202\estim_pen7_003_191004_144648\',...
               };
% Datafiles
ExplOD = {'CJ202\estim_pen7_003_191004_133556\CJ202_exp_datafile_023.mat',...
          'CJ202\estim_pen7_003_191004_144648\CJ202_exp_datafile_023.mat',...
         };
% Output directory
OutputlOD = {'CJ202\estim_pen7_004_conjoined\',...
             'CJ202_exp_datafile_023_conjoined.mat',...
            };
if isempty(dir([basePath OutputlOD{1}]))
    mkdir([basePath OutputlOD{1}]);
end
%% First, stitch the datafiles together
% Verify that the datafiles can be conjoined
Exp = cell(1,2); FN = cell(1,2);
for ei = 1:size(ExplOD,2)
    Exp{ei} = load([basePath ExplOD{ei}]);
    FN{ei} = fieldnames(Exp{ei});    
end
for ti = [1:8,10,11] % Magic number - the number of fields included in one of my experimental datafiles. There are 13 fields, but field 9 is the first cell
    if any(Exp{1}.(FN{1}{ti}) ~= Exp{2}.(FN{2}{ti}))
        disp('Incompatible datafiles');
        keyboard;
    end
end
for ti = [9,12,13] % Magic number - the number of fields included in one of my experimental datafiles. There are 13 fields, but field 9 is the first cell
    if any(cell2mat(Exp{1}.(FN{1}{ti})(2:end,1)) ~= cell2mat(Exp{2}.(FN{2}{ti})(2:end,1)))
        disp('Incompatible datafiles');
        keyboard;
    end
end
% All good - conjoin them!
Exp{3} = Exp{1};
Exp{3}.(FN{1}{5}) = Exp{3}.(FN{1}{5}) + Exp{2}.(FN{2}{5});
Exp{3}.(FN{1}{10}) = Exp{3}.(FN{1}{10}) + Exp{2}.(FN{2}{10});
Exp{3}.(FN{1}{11}) = [Exp{3}.(FN{1}{11}),Exp{2}.(FN{2}{11})];
Exp{3}.(FN{1}{12}) = [Exp{3}.(FN{1}{12});Exp{2}.(FN{2}{12})(2:end,:)];
Exp{3}.(FN{1}{13}) = [Exp{3}.(FN{1}{13});Exp{2}.(FN{2}{13})(2:end,:)];
% Save it out
AMP = Exp{3}.(FN{1}{1});
DUR = Exp{3}.(FN{1}{2});
E_DIAM = Exp{3}.(FN{1}{3});
E_MAT = Exp{3}.(FN{1}{4});
n_REP = Exp{3}.(FN{1}{5});
CHN = Exp{3}.(FN{1}{6});
PARAM = Exp{3}.(FN{1}{7});
PULSE = Exp{3}.(FN{1}{8});
E_MAP = Exp{3}.(FN{1}{9});
n_Trials = Exp{3}.(FN{1}{10});
rand_order = Exp{3}.(FN{1}{11});
StimParams = Exp{3}.(FN{1}{12});
TrialParams = Exp{3}.(FN{1}{13});
save([basePath OutputlOD{1} OutputlOD{2}],'AMP','DUR','E_DIAM','E_MAT','n_REP','CHN','PARAM','PULSE','E_MAP','n_Trials','rand_order','StimParams','TrialParams');
%% Time to stitch Intan together
% Processing will be done afterwards - this will allow egregious errors to
% be hopefully picked up by the processing scripts
amp1_fid = fopen([basePath IntanlOD{1} 'amplifier.dat'],'r');
amp2_fid = fopen([basePath IntanlOD{2} 'amplifier.dat'],'r');
dig1_fid = fopen([basePath IntanlOD{1} 'digitalin.dat'],'r');
dig2_fid = fopen([basePath IntanlOD{2} 'digitalin.dat'],'r');
out1_fid = fopen([basePath OutputlOD{1} 'amplifier.dat'],'W');
out2_fid = fopen([basePath OutputlOD{1} 'digitalin.dat'],'W');
nChn = 32; FS = 30000; T = 256; % Magic numbers.
% Amplifier waveforms
while ~feof(amp1_fid)
    f = fread(amp1_fid,[nChn FS*T],'int16');
    fwrite(out1_fid,f,'int16');
end
fclose(amp1_fid);
while ~feof(amp2_fid)
    f = fread(amp2_fid,[nChn FS*T],'int16');
    fwrite(out1_fid,f,'int16');
end
fclose(amp2_fid);
fclose(out1_fid);
% Digital Lines
while ~feof(dig1_fid)
    f = fread(dig1_fid,FS*T,'uint16');
    fwrite(out2_fid,f,'uint16');
end
fclose(dig1_fid);
while ~feof(dig2_fid)
    f = fread(dig2_fid,FS*T,'uint16');
    fwrite(out2_fid,f,'uint16');
end
fclose(dig2_fid);
fclose(out2_fid);
% Copy info.rhs as well
copyfile([basePath IntanlOD{1} 'info.rhs'],[basePath OutputlOD{1} 'info.rhs']);
%% Now for Blackrock
% Honestly, I haven't a clue how to do this while preserving the original
% Blackrock formats. Instead, this will operate on processed Blackrock
% files, and I'll just hope that nothing breaks.
denoise = dir([basePath BlackrocklOD{1} '*dnBR.dat']);
trig = dir([basePath BlackrocklOD{1} '*trigBR.dat']);
dn1_fid = fopen([basePath BlackrocklOD{1} denoise.name],'r');
tr1_fid = fopen([basePath BlackrocklOD{1} trig.name],'r');
denoise2 = dir([basePath BlackrocklOD{2} '*dnBR.dat']);
trig2 = dir([basePath BlackrocklOD{2} '*trigBR.dat']);
dn2_fid = fopen([basePath BlackrocklOD{2} denoise2.name],'r');
tr2_fid = fopen([basePath BlackrocklOD{2} trig2.name],'r');
out1_fid = fopen([basePath OutputlOD{1} denoise.name],'W');
out2_fid = fopen([basePath OutputlOD{1} trig.name],'W');
nChn = 96; FS = 30000; T = 256; % Magic numbers.
% Amplifier waveforms
while ~feof(dn1_fid)
    f = fread(dn1_fid,[nChn FS*T],'int16');
    fwrite(out1_fid,f,'int16');
end
fclose(dn1_fid);
while ~feof(dn2_fid)
    f = fread(dn2_fid,[nChn FS*T],'int16');
    fwrite(out1_fid,f,'int16');
end
fclose(dn2_fid);
fclose(out1_fid);
% Digital Lines
while ~feof(tr1_fid)
    f = fread(tr1_fid,FS*T,'double');
    fwrite(out2_fid,f,'double');
end
fclose(tr1_fid);
% At this point, we need to increment every number coming out of tr2 - 0 is
% actually the last timepoint in the previous file.
startTime = ((denoise.bytes/2)/nChn); % Start time in samples
while ~feof(tr2_fid)
    f = fread(tr2_fid,FS*T,'double');
    f = f + startTime;
    fwrite(out2_fid,f,'double');
end
fclose(tr2_fid);
fclose(out2_fid);
end