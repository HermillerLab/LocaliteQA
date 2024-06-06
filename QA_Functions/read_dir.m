function info = read_dir(dir_subj, subj, dir_QA)
% This function reads the directory of the subject 
% and extracts session information

% List all directories in the Sessions folder
timestamps=dir(fullfile(dir_subj,'Sessions/Session*')); 
timestamps={timestamps.name}; 

% Loop through each timestamp to extract session date and session name
for i=1:length(timestamps) 
    temp=timestamps{i}; 
    dates{i}=temp(9:16);
    try 
        % Try to parse the XML file to get session information
        xml=parseXML(fullfile(dir_subj,'Sessions/',timestamps{i},'Session.xml')); 
        sessions{i}=xml.Children(4).Children.Data;
    catch
        sessions{i}='empty_xml'; 
    end
end

% Combine session information into a matrix
info(:, 1:3) = [sessions', dates', timestamps']; 

% Save the session information matrix as 'info.mat' in the QA directory
save(fullfile(dir_QA,'info.mat'),'info');

end