function instrmt = read_instrmt(dir_subj, dir_QA, info, soi, tmstmp)
% This function reads the insturent marker data for each session

% Extract session names and timestamps from the 'info' matrix
sessions=info(:,1);  
timestamps=info(:,3);

% Loop through each session of interest (soi)
for d = 1:numel(soi)
    % Define the directory path where instrument marker data is located
    dir_data = fullfile(dir_subj, 'Sessions', timestamps{soi(d)}, 'InstrumentMarkers');

    % Open the most recent instrument marker XML file   
    files = dir(fullfile(dir_data, sprintf('InstrumentMarker*.xml'))); 
    xml = char({files(end).name});
    pxml = parseXML(fullfile(dir_data, xml));
    instrmt_timeset = xml(17:24);

    % Get instrument markers 4D matrix and save in rows
    ins = [];
    for i = 1:numel(pxml.Children)
        if strcmp(pxml.Children(i).Name, 'InstrumentMarker') == 1 && ...
           strcmp(pxml.Children(i).Attributes(3).Value, 'true') == 1
            ins(end+1, 1) = i;
        end
    end

    % Create a matrix to store instrument marker coordinates
    matrix = [];
    for i = 1:numel(ins)
        for m = 1:16
            matrix(i, m) = str2double(pxml.Children(ins(i)).Children(2).Children(2).Attributes(m).Value);
        end
    end
    
    % Store the instrument marker matrix under the session name
    instrmt.(sprintf('%s', sessions{soi(d)})) = matrix;
end

% Display instrument coordinates
disp('===========================================================================')
disp('                           INSTRUMENT COORDINATES                          ')
disp('===========================================================================')
disp(instrmt)

% Save the data and timeset in a .mat file in the QA_Results directory
save(fullfile(sprintf('%s/instmt_%s.mat', dir_QA, tmstmp)), 'instrmt','instrmt_timeset');

end