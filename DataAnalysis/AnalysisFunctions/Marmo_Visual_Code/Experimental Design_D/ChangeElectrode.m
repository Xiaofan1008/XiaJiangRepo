function ChangeElectrode()
    % Change electrode labels in StimParams and save new file
    E_Mapnumber_full = input('Please enter "6" for four shank, "5" for flexible single shank, "3" for rigid single shank: \n');
    
    % Prompt for file
    [DATAFILE, FOLDER] = uigetfile('*.mat', 'Select experimental .mat file');
    if isequal(DATAFILE, 0)
        disp('No file selected.');
        return;
    end
    cd(FOLDER);
    data = load(DATAFILE);
    
    % Check required variables
    if ~isfield(data, 'StimParams') || ~isfield(data, 'E_Mapnumber')
        error('Missing required variables in the data file.');
    end
    StimParams = data.StimParams;

    % Get electrode name mappings
    oldNames = input('Enter OLD channel names (comma-separated): ', 's');
    newNames = input('Enter NEW channel names (comma-separated): ', 's');
    oldNames = split(strtrim(oldNames), ',');
    newNames = split(strtrim(newNames), ',');
    
    if length(oldNames) ~= length(newNames)
        cerror('Number of old and new names must match.');
    end

    % Load electrode map
    E_MAP = ProbeMAP();
    E_MAP = E_MAP(:, E_Mapnumber_full); % Adjust based on E_Mapnumber
    if isempty(E_MAP{35})
        E_MAP=E_MAP(1:33);
    end

    % Helper function to map label to E_MAP name
    function mapped = mapElectrodeNames(nameList)
        mapped = cell(size(nameList));
        if E_Mapnumber_full == 6
            for i = 1:length(nameList)
                label = nameList{i};
                parts = regexp(label, 'P([A-D])S(\d)E(\d+)', 'tokens', 'once');
                if isempty(parts)
                    error('Invalid electrode name: %s', label);
                end
                port = parts{1};
                shank = str2double(parts{2});
                elec = str2double(parts{3});
                if isnan(shank) || isnan(elec)
                    error('Invalid shank or electrode in %s', label);
                end
                ch = elec;
                if shank == 3 || shank == 4
                    ch = ch + 16;
                elseif shank > 4
                    error('Shank number too large in %s', label);
                end
                switch port
                    case 'B', ch = ch + 32;
                    case 'C', ch = ch + 64;
                    case 'D', ch = ch + 96;
                end
                mapped{i} = E_MAP{ch + 1, 1}; % +1 for 1-based indexing
            end
        else
            for i = 1:length(nameList)
                ch = str2double(nameList{i});
                mapped{i} = E_MAP{ch + 1, 1};
            end
        end
    end

    % Map names
    mappedOld = mapElectrodeNames(oldNames);
    mappedNew = mapElectrodeNames(newNames);

    % Replace in StimParams
    for i = 2:size(StimParams, 1)
        for j = 1:length(mappedOld)
            if strcmp(StimParams{i,1}, mappedOld{j})
                StimParams{i,1} = mappedNew{j};
            end
        end
    end
    data.StimParams = StimParams;

    % Save to new file
    [~, baseName, ~] = fileparts(DATAFILE);
    newFileName = [baseName, '_changed_electrodes.mat'];
    save(newFileName, '-struct', 'data');
    disp(['Saved updated file as: ', newFileName]);
end


