 function stdStruct=stdstruct(stdStruct, ID, IDstruct,N)
check=['T', num2str(ID)];
stdStruct=[stdStruct, std(IDstruct.(check),[],2)/sqrt(N)];
 end