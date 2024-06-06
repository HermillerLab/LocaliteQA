clear all
addpath('QA_Functions');

% Display citation information
fprintf('PLEASE USE THE FOLLOWING CITATION WHEN REFERENCING OR UTILIZING THIS TOOLBOX\n')
fprintf('\nBurns, M. R., & Hermiller, M. S. (2024). LocaliteQA: https://github.com/HermillerLab/LocalityQA\n')
fprintf("\nTHANK YOU\n")

% Automatically detect the current path
current_path = pwd;
dir_subj = '';  % Initialize dir_subj to an empty string

% Check if folder contains valid files and folders for test
while isempty(dir_subj) || ~isfolder(dir_subj) || ...
      ~exist(fullfile(dir_subj, 'Sessions'), 'dir') || ...
      ~exist(fullfile(dir_subj, 'PatientData.xml'), 'file') || ...
      ~exist(fullfile(dir_subj, 'BinData'), 'dir')

    fprintf('\nSelect your participant folder (verify that selected folder contains a "Sessions" folder)\n');
    fprintf('If needed, see instructions for additional details on accessing correct subject folder\n');
    dir_subj = uigetdir(current_path, '');


    if isempty(dir_subj)
        fprintf('Folder selection canceled. Please try again.\n');
    elseif ~isfolder(dir_subj)
        fprintf('Invalid folder selected. Please select a valid participant folder.\n');
    elseif ~exist(fullfile(dir_subj, 'Sessions'), 'dir')
        fprintf('Selected folder does not contain a "Sessions" subfolder. Please select a different folder.\n');
    elseif ~exist(fullfile(dir_subj, 'PatientData.xml'), 'file')
        fprintf('Selected folder does not contain a "PatientData.xml" file. Please select a different folder.\n');
    elseif ~exist(fullfile(dir_subj, 'BinData'), 'dir')
        fprintf('Selected folder does not contain a "BinData" subfolder. Please select a different folder.\n');
    end
end

% Extract subject folder name 
slash_idx = strfind(dir_subj, filesep);
subj = dir_subj(slash_idx(end)+1:end); 
dir_root = dir_subj(1:slash_idx(end-1));

% Create directory for QA results
dir_arch = fullfile(current_path, 'QA_results');
no_break = false;

% Prompt the user for the folder name they want to create
prompt = 'Enter name for result folder: ';
user_folder_name = input(prompt, 's');


try 
    t = datetime('now'); 
    tmstmp = sprintf('%4g%02g%02g', t.Year, t.Month, t.Day); 
catch
    t = fix(datetime('now')); 
    tmstmp = sprintf('%4g%02g%02g', t(1), t(2), t(3)); 
end

folder_name = [user_folder_name '_' tmstmp];

% Create subject-specific directory for results
[status,msg] = mkdir(dir_arch, folder_name); 
dir_QA = fullfile(dir_arch, folder_name); 
info = read_dir(dir_subj, subj, dir_QA);


% Start logging session information
diary(fullfile(dir_QA, sprintf('QA_%s.txt', folder_name)));

% Display session information
disp('SESSIONS');
disp(info(:, 1)); 
disp('DATES'); 
disp(info(:, 2));

% Get the number of available sessions
num_sessions = size(info, 1);

% Prompt user to select sessions
valid_sessions = false;
while ~valid_sessions
    user_input = input('Enter session(s) (e.g., 1,2): ', 's');
    session_indices = strsplit(user_input, ',');
    session_indices = cellfun(@(x) str2double(strtrim(x)), session_indices);
    
    % Validate each session index
    valid_sessions = all(~isnan(session_indices) & session_indices >= 1 & session_indices <= num_sessions);
    
    if ~valid_sessions
        disp('Invalid input. Please enter valid session number(s).');
    end
end

% Convert session indices to numeric format
soi = session_indices;

% Parse XML files and gather information
[triggers, breaks] = read_triggers(dir_subj, subj, dir_QA, info, soi, tmstmp, no_break);
[entry, target] = read_entry(dir_subj, dir_QA, info, soi, tmstmp);
[instrmt] = read_instrmt(dir_subj, dir_QA, info, soi, tmstmp);
[RMSdev] = read_RMSdev(dir_subj, dir_QA, info, soi, tmstmp);

% Calculate deviations
[calcdev] = calc_dev(dir_subj, dir_QA, subj, info, soi, triggers, entry, instrmt, target, RMSdev, tmstmp, folder_name);

% Display citation information again
fprintf('PLEASE USE THE FOLLOWING CITATION WHEN REFERENCING OR UTILIZING THIS TOOLBOX\n')
fprintf('\nBurns, M. R., & Hermiller, M. S. (2024). LocaliteQA: https://github.com/HermillerLab/LocalityQA\n')
fprintf("\nTHANK YOU\n")

% Turn off diary logging
diary off


