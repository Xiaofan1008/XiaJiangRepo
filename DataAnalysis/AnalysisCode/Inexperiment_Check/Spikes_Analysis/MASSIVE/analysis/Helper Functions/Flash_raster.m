function count_sp_tr=Flash_raster(chn)
dbstop if error
%% Load in data from current directory

count_sp_tr=[];
plotHist=1;
plot_individualhist=1;
plot_hist_rasterstack=0;
plot_raster=1;
loadNREP;

fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir('info.rhs');
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=1;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end

try
    loadTrig(0);
catch
    cleanTrig_flash;%cleanTrig;%
end
d = Depth(E_Mapnumber); 
        
        %%%%%%%%%%%%%raster input
                chndat=zeros(nChn,401);
                trig=loadTrig(0);
                trig(trig==-500)=[];
                trig=trig(round(1:4.5:900));
                theseTrig = trig./30;
                nT=length(trig);
                %% Set up the raster data structure
                BIN = [-200 200];
                SPACING = 200; SMOOTHING = 2; MAX = 400;
                % Grab Data for the Inset Figure
                FS = 30000; file = pwd;
                X = BIN(1):1/30:BIN(2) - 1/30;
                mu = zeros(nT,(FS/1e3)*diff(BIN));
                for chn=chn(1):chn(end)
                    xdata = [];
                    ydata = [];
                    chn1 = d(chn);
                    sp = loadSpikes(chn1);
                for tr = 1:nT
                    theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                    count_sp_tr=[count_sp_tr length(theseSp)];
                    for i = 1:length(theseSp)
                        xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
                        ydata = [ydata, tr*(MAX/nT)];
                    end
                    % Also, extract MUA waveforms
                    %                 OFFSET = cast(nChn*2*(trig(tr)+(BIN(1)*FS/1e3)),'int64');
                    %                 fseek(mFID,OFFSET,'bof');
                    %                 v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
                    %                 mu(tr,:) = v(chn,:);
                    % Also, calculate phase for each trial
                    %     phase = generatePhaseVector(lfp,[4,5,6,10,20,80]);
                    %     for b = 1:6
                    %         thisPhase(tr,b) = phase{b}(cast(theseTrig(tr)/30 + offset(b),'int64'));
                    %     end
                end
                
                    Z = hist(xdata,0:400); %#ok<HIST>
                    chndat(chn1,:)=Z;
                end
                %% Plot the raster
                
                yScale = MAX./nT;
                if plot_hist_rasterstack==1
                    figure; hold on; axM = gca;
                    subplot(2,1,1)
                    hold on; axM = gca;
                    line([200 200],[0 MAX],'Color','b');
                    if plot_raster==1
                        for i = 1:size(xdata,2)
                            line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',1.25);
                        end
                    end
                else
                    figure; hold on; axM = gca;
%                     yline(MAX*size(relatedtrials_all,1)*counter_trials,'LineWidth',3);
%                     yline((counter_trials_row)*MAX+(counter_trials-1)*MAX*size(relatedtrials_all,1),'LineWidth',1);
%                     line([200 200],[0 MAX*size(relatedtrials_all,1)*size(relatedtrials_all,2)],'Color','b');
                     if plot_raster==1
                        for i = 1:size(xdata,2)
                            line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)],'Color','k','LineWidth',1.25);
                        end
                    end
                end
                
                Z = hist(xdata,0:400); %#ok<HIST>
            maxval=max(Z(150:250));
            if plotHist==1
                if plot_individualhist==1
                    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                    rate = (1000/nT)*conv(Z,window);
                    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                    plot(axM,rate,'b','LineWidth',2);
                else
                    sum_xhist=Z+sum_xhist;
                end
            end

xticks(axM,0:50:400);
xticklabels(axM,-200:50:200);
xline(200, 'r');
set(gca,'TickDir','out');
end