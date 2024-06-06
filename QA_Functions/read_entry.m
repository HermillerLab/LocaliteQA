function [entry, target] = read_entry(dir_subj, dir_QA, info, soi, tmstmp)
% This function reads the entry and target coordinates from XML files

% Extract session names and timestamps from the 'info' matrix
sessions=info(:,1); 
timestamps=info(:,3);

% Loop through each selected session of interest (soi)
for d = 1:numel(soi)
    % Construct the directory path for XML files
    dir_data = fullfile(dir_subj, 'Sessions', timestamps{soi(d)}, 'EntryTarget');

    % Find the latest XML file
    files = dir(fullfile(dir_data,sprintf('EntryTarget*.xml')));
    xml = char({files(end).name});
    pxml = parseXML(fullfile(dir_data, xml));
    
    % Initialize variable to store selected entry
    Selected = [];  
    % Find the entry node marked as 'true'
    for i = 1:numel(pxml.Children)
        if strcmp(pxml.Children(i).Name,'Entry')==1 && strcmp(pxml.Children(i).Attributes(3).Value,'true')==1
            Selected(end+1,1)=i;
        end    
    end
    
    %% Get entry coordinates
    ent = [];  
    for m = 1:3
        ent(1,m) = str2double(pxml.Children(Selected).Children(2).Children(2).Attributes(m).Value);
    end
    entry.(sprintf('%s',sessions{soi(d)})) = ent;
    
    %% Get target coordinates
    Tar_Selected=Selected+2;
    tar = [];  
    for m = 1:3
        tar(1,m) = str2double(pxml.Children(Tar_Selected).Children(2).Children(2).Attributes(m).Value);
    end
    target.(sprintf('%s',sessions{soi(d)})) = tar;
    
    
    %% Get rotation angle and reference
    Rot_Selected=Selected+4;
    rotation_ang.(sprintf('%s', sessions{soi(d)})) = str2double(pxml.Children(Rot_Selected).Attributes(1).Value);
    rot = [];  
    for m = 1:3
        rot(1,m) = str2double(pxml.Children(Rot_Selected).Children(2).Children(2).Attributes(m).Value);
    end
    rotation_ref.(sprintf('%s',sessions{soi(d)})) = rot;

end

% Display the extracted entry and target coordinates
disp('==========================================================================='); 
disp('                              ENTRY COORDINATES                            ');
disp('===========================================================================');
disp(entry)
disp('==========================================================================='); 
disp('                      PREDETERMINED TARGET COORDINATES                     ');
disp('===========================================================================');
disp(target)

end
