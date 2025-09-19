function ActiveLayerStruct=StructureLayergeneration
%used to load the .mat files and sum the responses based on the location of
%the stim electrode
D_data=dir;
fileptochange=pwd;
ActiveLayerStruct=[];
countnoexist=0;
numberelect4=0;
for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    cd(currD)
    try
        file=dir('ActivatedDepth*');
        load(file(1).name,'depthactivated')
    catch
        %doesn't exist
        countnoexist=countnoexist+1;
        cd(fileptochange)
        continue
    end
    
    for numsing=1:size(file,1)
        load(file(numsing).name,'depthactivated')
        for Elect=1:2
            layerstimone=find(cell2mat(depthactivated(2:end,Elect+1)));
            if ~isempty(layerstimone)
                
                if layerstimone==1 %layer 1
                    check=['L', num2str(layerstimone)];
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                elseif layerstimone==2 %layer 2/3
                    check='L2_3';
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                elseif layerstimone==3 %layer 4
                    check=['L', num2str(layerstimone+1)];
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                    
                elseif layerstimone==4 %layer 5
                    check=['L', num2str(layerstimone+1)];
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                    numberelect4=numberelect4+1;
                elseif layerstimone==5 %layer 6a
                    check='L6a';
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                    
                elseif layerstimone==6 %layer 6b
                    check='L6b';
                    if isfield(ActiveLayerStruct,check)
                        ActiveLayerStruct.(check) = ActiveLayerStruct.(check)+cell2mat(depthactivated(2:end,4));
                    else
                        ActiveLayerStruct.(check)= cell2mat(depthactivated(2:end,4));
                    end
                end
            end
        end
    end
    cd(fileptochange)
end
end

