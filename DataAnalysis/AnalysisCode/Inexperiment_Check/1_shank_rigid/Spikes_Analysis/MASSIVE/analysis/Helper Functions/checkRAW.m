function checkRAW
%% Variables
nChn = 32;
FS = 30000;
BIN = [-500 500; -946 946];
t = 1150;
%% Loading
loadTrig;
v_fid = fopen('amplifier.dat','r');
for b = 1:2
OFFSET = cast(nChn*2*(FS/1e3),'int64')*cast((trig(t)+BIN(b,1)),'int64');
fseek(v_fid,OFFSET,'bof');
v = fread(v_fid,[nChn, (FS/1e3)*diff(BIN(b,:))],'int16') .* 0.195;
%% Plotting
X = 1/30:1/30:size(v,2)/30;
figure; hold on
plot(X,v(14,:))
end
fclose(v_fid);
end
