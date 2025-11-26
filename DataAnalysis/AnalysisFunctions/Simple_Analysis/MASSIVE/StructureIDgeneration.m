 function IDstruct=StructureIDgeneration(values, ID, IDstruct)
 %creates Id Structure containing total number of spikes per trial for each
 %electrode grouped according to trial ID
 
 %INPUT - total number of spikes for this trial/repeat at each electrode 
 %(values), ID of the trial, ID structure to be concatenated
 
 %OUTPUT - ID STRUCTURE
 
    check=['T', num2str(ID)];
    if isfield(IDstruct,check)
        IDstruct.(check) = [IDstruct.(check), values];
    else
        IDstruct.(check)= values;
    end
 end