function stdStruct=stderrorstruct(stdStruct, ID, IDstruct)
check=['T', num2str(ID)];
if isfield(IDstruct, check)
    N=length(IDstruct.(check));
    stdStruct=[stdStruct, std(IDstruct.(check),[],2)/sqrt(N)];
else
    IDcell = struct2cell(IDstruct); %used to find number of electrodes
    stdStruct=[stdStruct, ones(size(IDcell{1},1),1)*(-500)];%if the requested trial does not exist then input impossible result to prevent any mixup and decrease number of required inputs
end

end