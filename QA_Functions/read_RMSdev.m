function rms_dev = read_RMSdev(dir_subj, dir_QA, info, soi, tmstmp)
% This function reads the RMS deviation data 
% at registration for each session

% Extract session information from 'info' matrix
sessions = info(:,1); 
dates = info(:,2); 
timestamps = info(:,3);

% Initialize flag to track if RMS deviations were found
rms_dev_found = false; 

for d = 1:numel(soi)
    % Define directory path for LocatorRegistrations
    dir_data = fullfile(dir_subj, 'Sessions', timestamps{soi(d)}, '/Registrations/LocatorRegistrations/');

    % Open the most recent Transformation XML file
    files_mat = dir(fullfile(dir_data, sprintf('Transformation%s*.xml', dates{soi(d)}))); 
    xml = char({files_mat(end).name});
    pxml = parseXML(fullfile(dir_data, xml));
    trans_timeset = xml(17:24);

    % Get transformation matrix and save in rows
    Transformation = [];
    for i = 1:numel(pxml.Children)
        if strcmp(pxml.Children(i).Name,'TransformationMatrix')
            Transformation(end+1,1) = i;
        end
    end

    matrix = [];
    for i = 1:numel(Transformation)
        for m = 1:16
            matrix(i, m) = str2num(pxml.Children(Transformation(i)).Children(2).Attributes(m).Value);
        end
    end

    trans_mat.(sprintf('%s', sessions{soi(d)})) = matrix;

    % open Transformation.txt file
    try 
        files = dir(fullfile(dir_data, 'Transformation*.txt'));
        last = char({files(end).name});
        txt = dir(fullfile(dir_data, last));
        filename = char(txt.name);

        file = fopen(fullfile(dir_data, filename), 'r');
        tline = fgetl(file); 
        rms = [];
        rot = [];
        transl = [];

        % Read data from the Transformation.txt file
        while ischar(tline)
            if ~isempty(strfind(tline, 'RMS Deviation'))
                rms = [rms regexp(tline, '\d*\.\d*', 'match')];
            end
            if ~isempty(strfind(tline, 'Rotation'))
                rot = [rot regexp(tline, '\d*\.\d*', 'match')];
            end
            if ~isempty(strfind(tline,'Translation'))
                transl = [transl regexp(tline,'\d*\.\d*','match')];
            end

            tline = fgetl(file);
        end

        fclose(file);

        if ~isempty(rms)
            % RMS deviations found
            rms_dev_found = true;
            rms_dev.(sprintf('%s', sessions{soi(d)})).landmark = rms(1);
            rms_dev.(sprintf('%s', sessions{soi(d)})).surface = rms(2);
            rms_dev.(sprintf('%s', sessions{soi(d)})).rotation = rot(1);
            rms_dev.(sprintf('%s', sessions{soi(d)})).translation = transl(1);
            rmsdev(d, :) = rms;
        else
            % No RMS deviations found
            disp(['No RMS deviations found for session: ' sessions{soi(d)}]);
        end

    catch
        rms_dev.(sprintf('%s', sessions{soi(d)})).landmark = NaN;
        rms_dev.(sprintf('%s', sessions{soi(d)})).surface = NaN;
    end
end

% Display RMS deviations at registration
disp('==========================================================================='); 
disp('                       RMS DEVIATIONS AT REGISTRATION                      ');
disp('===========================================================================');

% Print message if no RMS deviations were found
if rms_dev_found
    disp(rmsdev');
else
    disp('No RMS deviations found.');
end

% Save RMS deviation data and transformation matrix in a .mat file in the QA directory
save(fullfile(sprintf('%s/RMSdv_%s.mat', dir_QA, tmstmp)), 'rms_dev', 'trans_mat');

end

