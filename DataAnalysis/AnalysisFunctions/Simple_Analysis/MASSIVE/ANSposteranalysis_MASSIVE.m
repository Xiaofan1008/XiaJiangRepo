%%ANS poster analysis MASSIVE

%% 1. Raster and PSTH
%loop folders

D_data=dir;
respsaveV1=[];
ratesaveV1=[];
xdatasaveV1=[];
ydatasaveV1=[];
VAsigcountV1=0;
respsaveV2=[];
ratesaveV2=[];
xdatasaveV2=[];
ydatasaveV2=[];
VAsigcountV2=0;
spkarraysaveV2=[];
spkarraysaveV1=[];
chn_range=[1 128];



for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
    catch
        continue
    end
    amplifier_channels=read_Intan_RHS2000_file;
    nChn=length(amplifier_channels);
    filepath=pwd;
    [filepathm,name,ext] = fileparts(filepath);
    ti=loadTrialInfo;
    trial=cell2mat(ti([false; cell2mat(ti(2:end,18))==10],1));
    load([name(1:end-14) '.sp.mat'],'sp')

    order=Depth(1);
    for ID=trial'%5:5:size(ti,1)-1
        
        for chn=chn_range(1):chn_range(end) %=1:64%65:128
            [resp,rate,xdata,ydata]=ANS_raster(ID,sp{chn});
            if ~isempty(resp) && any(resp>100) && chn<65
                %title(['chn: ' num2str(find(chn==order)) ' ID: ' num2str(ID)])
                %saveas(gcf,['Plot' num2str(chn) '_' num2str(ID) '.fig']);
                %respsaveV2=[respsaveV2 resp];
               
               [~,I]= max(rate);
               timestart=I-100-5;
               if timestart<2
                   timeend=12;
                   timestart=2;
               else
                   timeend=I-100+5;
               end
                [spkarrayV2]=ANS_SigmoidGenerator(chn,nChn,timestart,timeend,ID,sp);
                                  halfmaxV2=max(spkarrayV2,[],2)/2;
                [~,V2_I]=min(abs(spkarrayV2-halfmaxV2));
                if V2_I==1
                    continue
                end
                 ratesaveV2=[ratesaveV2; rate];
                  spkarraysaveV2=[spkarraysaveV2;spkarrayV2];
                %xdatasaveV2=[xdatasaveV2,xdata];
                %ydatasaveV2=[ydatasaveV2,ydata];
                VAsigcountV2=VAsigcountV2+1;
            elseif ~isempty(resp) && any(resp>100) && chn>64
                %respsaveV1=[respsaveV1 resp];
                [~,I]= max(rate);
                timestart=I-100-5;
                if timestart<2
                    timeend=12;
                    timestart=2;
                else
                    timeend=I-100+5;
                end
                [spkarrayV1]=ANS_SigmoidGenerator(chn,nChn,timestart,timeend,ID,sp);
                halfmaxV1=max(spkarrayV1,[],2)/2;
                [~,V1_I]=min(abs(spkarrayV1-halfmaxV1));
                if V1_I==1
                    continue
                end
                 spkarraysaveV1=[spkarraysaveV1;spkarrayV1];
                  ratesaveV1=[ratesaveV1; rate];
                %xdatasaveV1=[xdatasaveV1,xdata];
                %ydatasaveV1=[ydatasaveV1,ydata];
                VAsigcountV1=VAsigcountV1+1;
            end
        end

        
    end
        fprintf(currD)
        fprintf(' completed\n')
        cd([D_data(k).folder])
        save('summarydat.mat','ratesaveV1','ratesaveV2','spkarraysaveV1','spkarraysaveV2')
    
end
