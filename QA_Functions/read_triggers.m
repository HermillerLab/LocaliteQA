function [triggers, breaks0] = read_triggers(dir_subj, subj, dir_QA, info, soi, tmstmp, no_break)
% This function reads trigger markers from Localite TMS session

% Extract session names, dates, and timestamps from the 'info' matrix
sessions = info(:, 1); 
dates = info(:, 2); 
timestamps = info(:, 3);

session_break_limit = 20; % minimum change to elicit a session break

disp('============================================================================'); 
disp('> > > >  SELECT SESSION TRIGGER FILE (LIKELY THE LARGEST LAST FILE)  < < < <'); 
disp('============================================================================'); 

for d = 1:numel(soi)
    dir_data = fullfile(dir_subj, 'Sessions', timestamps{soi(d)}, 'TMSTrigger');
    files = dir(fullfile(dir_data, sprintf('*Coil1_%s*.xml', dates{soi(d)}))); 
    file_size = {files.bytes}; 
    fprintf('%s: \n', sessions{soi(d)}); 
    disp(file_size);
    reply = str2num(input('Press ENTER to use last file in session (or type file number "{[####]}"): ','s')); 

    if isempty(reply) 
        reply = numel(files); 
    end
    use{d} = reply;
end

disp('parsing trigger xml files');

%% open TriggerMarkers_Coil1_timestamp.xml

breaks = 0;

for d = 1:numel(soi)
    invalid_trig_file = true;
    num_try = 1;
    while invalid_trig_file
        tic
        dir_data = fullfile(dir_subj, 'Sessions', timestamps{soi(d)}, 'TMSTrigger');
        files = dir(fullfile(dir_data, sprintf('*Coil0_%s*.xml', dates{soi(d)})));
        x = use{d};
        coil0_xml = char({files(x).name});
        coil1_xml = strrep(coil0_xml,'Coil0','Coil1');
        coil0 = parseXML(fullfile(dir_data, coil0_xml));
        try
            coil1 = parseXML(fullfile(dir_data, coil1_xml));
        catch
            fprintf(2, '\nUnable to process coil1 datafile: %s\n Please check that correct trigger file was selected\n\n', coil1_xml)
            coil1 = coil0;
        end

        % Get trigger markers 4D matrix and save in rows
        TriggerMarker0 = [];
        for i = 1:numel(coil0.Children)
            if strcmp(coil0.Children(i).Name, 'TriggerMarker') == 1
                TriggerMarker0(end+1, 1) = i;
            end
        end
        TriggerMarker1 = [];
        for i = 1:numel(coil1.Children)
            if strcmp(coil1.Children(i).Name, 'TriggerMarker') == 1
                TriggerMarker1(end+1, 1) = i;
            end
        end

              if numel(TriggerMarker0) == 0 && numel(TriggerMarker1) == 0
            num_try = num_try + 1;
            if num_try < 100
                fprintf(2, 'No trigger markers found. Attempts left: %d\n', 100 - num_try);
                fprintf(2, 'The selected trigger file for session %s is empty.\nPlease choose another trigger file or hit ctrl-c to exit\n\n',...
                    sessions{soi(d)});
                file_size = {files.bytes}; 
                fprintf('%s:', sessions{soi(d)}); 
                disp(file_size);
                reply = str2num(input('Press ENTER to use last file in session (or type file number "{[####]}"): ','s'));
                if isempty(reply) 
                    reply = numel(files); 
                end
                use{d} = reply;
            else
                fprintf(2, 'Final attempt unsuccessful. Exiting program.\n\n');
                exit;
            end
        else
            invalid_trig_file = false;
        end

    end

     % di/dt (A/us) for each trigger
    didt0 = [];
    part_session = 1;
    breaks0 = [0];
    ps = 0;
    for i = 1:numel(TriggerMarker0)
        ps = ps + 1;
        didt0(i) = str2num(coil0.Children(TriggerMarker0(i)).Children(2).Children(2).Attributes(2).Value);

        if ps > 20 && ~no_break
            ma_didt = nanmean(didt0(i-20:i-1)); % moving average of the last 20 items
            if abs(didt0(i) - ma_didt) > session_break_limit
                part_session = part_session + 1;
                ps = 0;
                breaks0(1) = breaks0(1) + 1;
                breaks0(breaks0(1)+1) = i;
            end
        end
    end

    didt1 = [];
    breaks1 = [0];
    part_session = 1;
    ps = 0;
    for i = 1:numel(TriggerMarker1)
        ps = ps + 1;
        didt1(i) = str2num(coil1.Children(TriggerMarker1(i)).Children(2).Children(2).Attributes(2).Value);

        if ps > 20 && ~no_break
            ma_didt1 = nanmean(didt1(i-20:i-1)); % moving average of the last 20 items
            if abs(didt1(i) - ma_didt1) > session_break_limit
                part_session = part_session + 1;
                ps = 0;
                breaks1(1) = breaks1(1) + 1;
                breaks1(breaks1(1)+1) = i;
            end
        end
    end
    
    % 4Dmatrix of each trigger
    matrix0 = [];
    for i = 1:numel(TriggerMarker0)
        for m = 1:16
            matrix0(i,m) = str2double(coil0.Children(TriggerMarker0(i)).Children(4).Attributes(m).Value); 
            if matrix0(i,m) == 0 
                matrix0(i,m) = NaN; 
            end 
        end
    end
    matrix1 = [];
    for i = 1:numel(TriggerMarker1)
        for m = 1:16
            matrix1(i,m) = str2double(coil1.Children(TriggerMarker1(i)).Children(4).Attributes(m).Value); 
            if matrix1(i,m) == 0 
                matrix1(i,m) = NaN; 
            end 
        end
    end

    % Record time of triggers
    time = [];
    for i = 1:numel(TriggerMarker0)
        time(i) = str2num(coil0.Children(TriggerMarker0(i)).Attributes(3).Value);
    end

    % Record amplitude for each trigger
    amp0 = []; 
    for i = 1:numel(TriggerMarker0)
        amp0(i) = str2num(coil0.Children (TriggerMarker0(i)).Children(2).Children(6).Attributes(2).Value);
    end
    amp1 = []; 
    for i = 1:numel(TriggerMarker1)
        amp1(i) = str2num(coil1.Children(TriggerMarker1(i)).Children(2).Children(6).Attributes(2).Value);
    end


    % Determine which coil file to use for each measure

    if d == 1
        disp('==========================================================================='); 
        disp('                              COIL FILES USED                              ');
        disp('===========================================================================');
    end
    if nanmean(didt1) > nanmean(didt0)

        didt = didt1;
        if length(breaks1) > size(breaks, 2)
            breaks = padarray(breaks, [0, length(breaks1) - size(breaks, 2)], NaN, 'post');
        end
        breaks(d, 1:length(breaks1)) = breaks1;
        fprintf('Coil0 used for di/dt for session %s\n', sessions{soi(d)})
    else
        didt = didt0;
        if length(breaks0) > size(breaks, 2)
            breaks = padarray(breaks, [0, length(breaks0) - size(breaks, 2)], NaN, 'post');
        end
        breaks(d, 1:length(breaks0)) = breaks0;
        if isnan(nanmean(didt0)) || nanmean(didt0) == 0
            fprintf(2, 'No valid didts were found so Coil1 was used for session %s\n',...
                sessions{soi(d)})
        else
            fprintf('Coil1 used for di/dt for session %s\n', sessions{soi(d)})
        end
    end

    if nansum(abs(nanmean(matrix0))) <= 4
        if nansum(abs(nanmean(matrix1))) <= 4
            matrix = matrix0;
            fprintf(2, 'No valid locations were found so Coil1 was used for session %s\n',...
                sessions{soi(d)})
        else
            matrix = matrix1;
            fprintf('Coil1 used for location matrix for session %s\n', sessions{soi(d)})
        end
    else
        matrix = matrix0;
        fprintf('Coil0 used for location matrix for session %s\n', sessions{soi(d)})
    end

    if any(isnan(amp0)) || all(amp0 == 0)
        amp = amp1;
        if isempty(amp1)
            fprintf(2, 'No valid amplitudes were found so Coil1 was used for session %s\n',...
                sessions{soi(d)});
        else
            fprintf('Coil1 used for amp for session %s\n', sessions{soi(d)});
        end
    else
        if nanmean(amp1) > nanmean(amp0)
            amp = amp1;
            fprintf('Coil1 used for amp for session %s\n', sessions{soi(d)});
        else
            amp = amp0;
            fprintf('Coil0 used for amp for session %s\n', sessions{soi(d)});
        end
    end
    
    if d == numel(soi)
        fprintf('\n')
    end
   
    % Break up and store the values
    if breaks(d, 1) == 0
        di_dt.(sprintf('%s', sessions{soi(d)})) = didt;
        ave_didt(:, d, 1) = nanmean(didt);
        triggers.(sprintf('%s', sessions{soi(d)})) = matrix; 
        ave_pulse(:, d, 1) = nanmean(matrix)';
        trig_times.(sprintf('%s', sessions{soi(d)}))=time;
        amplitude.(sprintf('%s',sessions{soi(d)})) = amp;
        valid_input(d, 1) = 1;

    elseif breaks(d, 1) > 0
        for i = 1:(breaks(d, 1)+1)
            if i == 1
                di_dt.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = didt(1: breaks(d, i+1)-1);
                ave_didt(:, d, i) = nanmean(didt(1: breaks(d, i+1)-1));
                triggers.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = matrix(1: breaks(d, i+1)-1, :); 
                ave_pulse(:, d, i) = nanmean(matrix(1: breaks(d, i+1)-1, :));
                trig_times.(sprintf('%s_%.0f', sessions{soi(d)}, i))=time(1: breaks(d, i+1)-1);
                amplitude.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = amp(1: breaks(d, i+1)-1);
            elseif i == breaks(d, 1)+1
                di_dt.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = didt(breaks(d, i):end);
                ave_didt(:, d, i) = nanmean(didt(breaks(d, i):end));
                triggers.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = matrix(breaks(d, i):end, :); 
                ave_pulse(:, d, i) = nanmean(matrix(breaks(d, i):end, :));
                trig_times.(sprintf('%s_%.0f', sessions{soi(d)}, i))=time(breaks(d, i):end);
                amplitude.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = amp(breaks(d, i):end);
            else
                di_dt.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = didt(breaks(d, i): breaks(d, i+1)-1);
                ave_didt(:, d, i) = nanmean(didt(breaks(d, i): breaks(d, i+1)-1));
                triggers.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = matrix(breaks(d, i): breaks(d, i+1)-1, :); 
                ave_pulse(:, d, i) = nanmean(matrix(breaks(d, i): breaks(d, i+1)-1, :))';
                trig_times.(sprintf('%s_%.0f', sessions{soi(d)}, i))=time(breaks(d, i): breaks(d, i+1)-1);
                amplitude.(sprintf('%s_%.0f', sessions{soi(d)}, i)) = amp(breaks(d, i): breaks(d, i+1)-1);
            end

            valid_input(d, i) = 1;
        end
    end

end

% Organize information for output
c = 1;
for i = 1:size(ave_pulse, 2)
    for j = 1:size(ave_pulse, 3)
        if valid_input(i, j) ~= 0
            ave_pulse_4disp(c, :) = ave_pulse(:, i, j)';
            c = c+1;
        end
    end
end

c = 1;
for i = 1:size(ave_didt, 2)
    for j = 1:size(ave_didt, 3)
        if valid_input(i, j) ~= 0
            ave_didt_4disp(c, :) = ave_didt(:, i, j)';
            c = c+1;
        end
    end
end

% Display averages 
disp('==========================================================================='); 
disp('                          AVERAGE TRIGGER LOCATION                         ');
disp('==========================================================================='); 
disp(rem_NaN_rotmax(ave_pulse_4disp));
disp('==========================================================================='); 
disp('                           AVERAGE RECORDED di/dt                          ');
disp('===========================================================================');
disp(ave_didt_4disp);

% Save .mat files 
save(fullfile(sprintf('%s/pulse_%s.mat', dir_QA, tmstmp)), 'triggers', 'trig_times', 'ave_pulse');
save(fullfile(sprintf('%s/didt_%s.mat', dir_QA, tmstmp)), 'di_dt', 'amplitude', 'ave_didt');

end

