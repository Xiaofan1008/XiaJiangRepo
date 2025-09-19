function avgStruct=meanstruct(avgStruct, ID, IDstruct)
check=['T', num2str(ID)];
if isfield(IDstruct, check)
    avgStruct=[avgStruct, mean(IDstruct.(check),2)];
else
    IDcell = struct2cell(IDstruct);%used to find number of electrodes
    avgStruct=[avgStruct, ones(size(IDcell{1},1),1)*(-500)];%if the requested trial does not exist then input impossible result to prevent any mixup and decrease number of required inputs
end

end