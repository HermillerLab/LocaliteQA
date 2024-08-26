function calcdev = calc_dev(dir_subj, dir_QA, subj, info, soi, triggers, entry, instrmt, target, RMSdev, tmstmp, folder_name,N,T)


% Extract session names and trigger/instrument fields
sessions = info(:,1);
trigger_fields = fieldnames(triggers);
inst_fields = fieldnames(instrmt);
calcdev = struct;
to_graph = 1:length(trigger_fields);

% Loop through each trigger field
for t = 1:length(trigger_fields)
    % Load trigger data
    trg{t} = triggers.(trigger_fields{t});
    num_trig = size(trg{t}, 1);
    % Clean trigger file for problematic triggers
    t2r = 0;
    for i = 1:num_trig
       if any(isnan(trg{t}(i,1:12)))
           t2r(1) = t2r(1) + 1;
           t2r(end+1) = i;
       end
    end
    
    % Remove triggers with NaN values
    if t2r(1) > 0        
        warning off backtrace
        warning('\n%.0f triggers removed for session %s due to NaNs, please look at trigger files\n', t2r(1), trigger_fields{t})
        warning on backtrace
        trg{t}(t2r(2:end),:) = [];  
    end
    
    % Handle case where all triggers are removed
    if t2r(1) == num_trig 
        fprintf(2, '\nAll triggers removed for session %s, please make sure correct session was selected and re-run\n\n', trigger_fields{t})
        to_graph(t) = NaN;
        calcdev.max_trg2prevtrg(t, :) = NaN(16,1);
        calcdev.ave_trg2prevtrg(t, :) = NaN(16,1);
        calcdev.var_pls2prevpls(t, :) = NaN(16,1);
        var_trg(t, :) = NaN(16,1);
        var_pls(t, :) = NaN(3,1);
        var_pls_dist0(t, 1) = NaN;
        imk2ent(t, 1) = NaN;
        ave_pls2ent(t, 1) = NaN;
        std_pls2ent(t, 1) = NaN;
        var_pls2ent(t, 1) = NaN;
        ave_pls2imk(t, 1) = NaN;
        std_pls2imk(t, 1) = NaN;
        var_pls2imk(t, 1) = NaN;
        
    else    
        % Find matching session for trigger session
        for i = 1:length(inst_fields)
            if strcmp(inst_fields{i}, trigger_fields{t}(1:length(inst_fields{i})))
                s = i;
                break
            elseif i == length(inst_fields)
                fprintf(2, 'ERROR: Cannot find match for session file %s\n', ...
                    trigger_fields{t})
            end
        end
        ins(t, :) = instrmt.(sprintf('%s', sessions{soi(s)}));
        ent(t, :) = entry.(sprintf('%s', sessions{soi(s)}));
        tar(t, :) = target.(sprintf('%s', sessions{soi(s)}));
        ent2tar(t,:) = tar(t,:)-ent(t,:);

        % Convert 4D trigger matrix to 3D vector
        for i = 1:numel(trg{t}(:, 1))
            pls{t}(i, 1) = trg{t}(i, 4); 
            pls{t}(i, 2) = trg{t}(i, 8); 
            pls{t}(i, 3) = trg{t}(i, 12); 
            pls{t}(i, 4) = trg{t}(i, 16); 
        end

        % Convert 4D instrument matrix to 3D vector
        % Get position info for instrument marker
        imk(t, 1) = ins(t, 4); 
        imk(t, 2) = ins(t, 8); 
        imk(t, 3) = ins(t, 12); 
        imk(t, 4) = ins(t, 16);

        % Get rotational info for instrument marker
        ins_rot{t}(1,1)=ins(t, 1);
        ins_rot{t}(1,2)=ins(t, 2);
        ins_rot{t}(1,3)=ins(t, 3);
        ins_rot{t}(1,4)=ins(t, 4);
        ins_rot{t}(2,1)=ins(t, 5);
        ins_rot{t}(2,2)=ins(t, 6);
        ins_rot{t}(2,3)=ins(t, 7);
        ins_rot{t}(2,4)=ins(t, 8);
        ins_rot{t}(3,1)=ins(t, 9);
        ins_rot{t}(3,2)=ins(t, 10);
        ins_rot{t}(3,3)=ins(t, 11);
        ins_rot{t}(3,4)=ins(t, 12);
        ins_rot{t}(4,1)=ins(t, 13);
        ins_rot{t}(4,2)=ins(t, 14);
        ins_rot{t}(4,3)=ins(t, 15);
        ins_rot{t}(4,4)=ins(t, 16);

        % Calculate deviations and variance (Euclidean distance; position changes only

        % Euclidean Distance between instrument and entry
        imk2ent(t, :) = sqrt((imk(t, 1) - ent(t, 1))^2 + (imk(t, 2) - ent(t,2))^2 + (imk(t,3) - ent(t,3))^2);
  
        % Euclidean Distance between instrument and target
        imk2tar(t, :) = sqrt((imk(t, 1) - tar(t, 1))^2 + (imk(t, 2) - tar(t,2))^2 + (imk(t,3) - tar(t,3))^2);
        
        % Euclidean Distance between target and entry
        tar2ent(t, :) = sqrt((tar(t, 1) - ent(t, 1))^2 + (tar(t, 2) - ent(t,2))^2 + (tar(t,3) - ent(t,3))^2);
        
        
        for p = 1:length(pls{t}(:, 1))
            % Euclidean Distance between pulses and entry, instrument, and  target
            pls2ent(t, p) = sqrt(((pls{t}(p, 1) - ent(t, 1))^2) + ((pls{t}(p, 2) - ent(t, 2))^2) + ((pls{t}(p, 3) - ent(t, 3))^2));
            pls2imk(t, p) = sqrt((pls{t}(p, 1) - imk(t, 1))^2 + (pls{t}(p, 2) - imk(t, 2))^2 + (pls{t}(p, 3) - imk(t, 3))^2);           
            pls2tar(t, p) = sqrt((pls{t}(p, 1) - tar(t, 1))^2 + (pls{t}(p, 2) - tar(t, 2))^2 + (pls{t}(p, 3) - tar(t, 3))^2);

            % Get rotational info for each trigger
            trg_rot{t}{p}(1,1)=trg{t}(p, 1);
            trg_rot{t}{p}(1,2)=trg{t}(p, 2);
            trg_rot{t}{p}(1,3)=trg{t}(p, 3);
            trg_rot{t}{p}(1,4)=trg{t}(p, 4);
            trg_rot{t}{p}(2,1)=trg{t}(p, 5);
            trg_rot{t}{p}(2,2)=trg{t}(p, 6);
            trg_rot{t}{p}(2,3)=trg{t}(p, 7);
            trg_rot{t}{p}(2,4)=trg{t}(p, 8);
            trg_rot{t}{p}(3,1)=trg{t}(p, 9);
            trg_rot{t}{p}(3,2)=trg{t}(p, 10);
            trg_rot{t}{p}(3,3)=trg{t}(p, 11);
            trg_rot{t}{p}(3,4)=trg{t}(p, 12);
            trg_rot{t}{p}(4,1)=0;
            trg_rot{t}{p}(4,2)=0;
            trg_rot{t}{p}(4,3)=0;
            trg_rot{t}{p}(4,4)=trg{t}(p, 16);
            
            % Use the Procrustes transformation returned by procrustes to analyze how it superimposes the triggers onto the instrument marker
            [d,z,transform]=procrustes(ins_rot{t},trg_rot{t}{p},"scaling",false,"reflection",false);
            procrustes_trg2ins{t}(p,1)=d; %dispartity between two matrices
            
            % Calculate the Euler angles from the rotation matrix
            R = transform.T;
            
            % Check for gimbal lock (is sy close to zero?); if so, then ignore yaw
            sy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1));
            singular = sy < 1e-6;
            
            if ~singular
                pitch{t}(p) = atan2d(R(3,2), R(3,3)); % tilting the handle forward/backward
                roll{t}(p) = atan2d(-R(3,1), sy); % rotating the coil around the tail (circling around the hotspot; aligning the "T" of the instrument marker)
                yaw{t}(p) = atan2d(R(2,1), R(1,1)); % side-to-side (butterfly arms flapping; see-saw)
            else
                pitch{t}(p) = atan2d(-R(2,3), R(2,2));
                roll{t}(p) = atan2d(-R(3,1), sy);
                yaw{t}(p) = 0;
            end
            
            % Calculate "target" of each pulse
            tar_each_pls{t}(p, 1)= trg_rot{t}{p}(1, 4) + pls2tar(t, p) * trg_rot{t}{p}(1, 1);
            tar_each_pls{t}(p, 2)= trg_rot{t}{p}(2, 4) + pls2tar(t, p) * trg_rot{t}{p}(2, 1);
            tar_each_pls{t}(p, 3)= trg_rot{t}{p}(3, 4) + pls2tar(t, p) * trg_rot{t}{p}(3, 1);
            
            % Distance each "target" is away from intended target
            dist_pls_tar_off_tar(t, p)= sqrt(((tar_each_pls{t}(p, 1) - tar(t, 1))^2) + ((tar_each_pls{t}(p, 2) - tar(t, 2))^2) + ((tar_each_pls{t}(p, 3) - tar(t, 3))^2));
                       
        end

        ave_pls_tar(t, :) = mean(tar_each_pls{t});

        % Count pulses N+ off entry and instrument and T+ off target
        pls2ent_OFF = (pls2ent(t, :) >= N); pls2ent_ctOFF(t, 1) = sum(pls2ent_OFF);
        pls2imk_OFF = (pls2imk(t,:) >= N); pls2imk_ctOFF(t, 1) = sum(pls2imk_OFF);
        dist_pls_tar_off_tar_OFF = (dist_pls_tar_off_tar(t,:) >= T); dist_pls_tar_off_tar_ctOFF(t, 1) = sum(dist_pls_tar_off_tar_OFF);     
    
        % Summary of pulses from entry, instrument, and target
        ave_pls2ent(t, 1) = mean(pls2ent(t, 1:length(pls{t}(:, 1))));
        std_pls2ent(t, 1) = std(pls2ent(t, 1:length(pls{t}(:, 1))));
        var_pls2ent(t, 1) = var(pls2ent(t, 1:length(pls{t}(:, 1))));
        ave_pls2imk(t, 1) = mean(pls2imk(t, 1:length(pls{t}(:, 1))));
        std_pls2imk(t, 1) = std(pls2imk(t, 1:length(pls{t}(:, 1))));
        var_pls2imk(t, 1) = var(pls2imk(t, 1:length(pls{t}(:, 1))));        
        ave_pls2tar(t, 1) = mean(pls2tar(t, 1:length(pls{t}(:, 1))));
        std_pls2tar(t, 1) = std(pls2tar(t, 1:length(pls{t}(:, 1))));
        var_pls2tar(t, 1) = var(pls2tar(t, 1:length(pls{t}(:, 1))));
        ave_dist_pls_tar_off_tar(t, 1) = mean(dist_pls_tar_off_tar(t, 1:length(pls{t}(:, 1))));
        std_dist_pls_tar_off_tar(t, 1) = std(dist_pls_tar_off_tar(t, 1:length(pls{t}(:, 1))));
        var_dist_pls_tar_off_tar(t, 1) = var(dist_pls_tar_off_tar(t, 1:length(pls{t}(:, 1))));
        ave_procrustes_trg2ins(t,1) = mean(procrustes_trg2ins{t}(:,1));
        std_procrustes_trg2ins(t,1) = std(procrustes_trg2ins{t}(:,1));
        var_procrustes_trg2ins(t,1) = var(procrustes_trg2ins{t}(:,1));
        ave_pitch(t,1) = mean(pitch{t}(:));
        std_pitch(t,1) = std(pitch{t}(:));
        var_pitch(t,1) = var(pitch{t}(:));
        ave_roll(t,1) = mean(roll{t}(:));
        std_roll(t,1) = std(roll{t}(:));
        var_roll(t,1) = var(roll{t}(:));
        ave_yaw(t,1) = mean(yaw{t}(:));
        std_yaw(t,1) = std(yaw{t}(:));
        var_yaw(t,1) = var(yaw{t}(:));
        
        % Displacement from previous
        for u = 1:p-1 
            trg2prevtrg{t}(u, :) = (trg{t}((u+1), :)) - (trg{t}(u, :));
            pls2prevpls(t, u) = sqrt((pls{t}(u+1, 1) - pls{t}(u, 1))^2 + (pls{t}(u+1, 2) - pls{t}(u, 2))^2 + (pls{t}(u+1, 3) - pls{t}(u, 3))^2);

            [d,z,transform]=procrustes(trg_rot{t}{u+1},trg_rot{t}{u},"scaling",false,"reflection",false);
            procrustes_trgvar{t}(u,1)=d; % dispartity between two matrice
            
            % Calculate the Euler angles from the rotation matrix
            R = transform.T;
            
            % Check for gimbal lock (is sy close to zero?); if so, then ignore yaw
            sy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1));
            singular = sy < 1e-6;
            
            if ~singular
                pitch_trgvar(t,u) = atan2d(R(3,2), R(3,3)); % tilting the handle forward/backward
                roll_trgvar(t,u) = atan2d(-R(3,1), sy); % rotating the coil around the tail (circling around the hotspot; aligning the "T" of the instrument marker)
                yaw_trgvar(t,u) = atan2d(R(2,1), R(1,1)); % side-to-side (butterfly arms flapping; see-saw)
            else
                pitch_trgvar(t,u) = atan2d(-R(2,3), R(2,2));
                roll_trgvar(t,u) = atan2d(-R(3,1), sy);
                yaw_trgvar(t,u) = 0;
            end
        end
         
        var_trg(t, :) = var(trg{t}); 
        var_pls(t, :) = var(pls{t}); 
        var_pls_dist0(t, 1) = sqrt(var_pls(t, 1)^2 + var_pls(t, 2)^2 + var_pls(t, 3)^2);

        calcdev.max_trg2prevtrg(t, :) = max(trg2prevtrg{t});
        calcdev.ave_trg2prevtrg(t, :) = mean(trg2prevtrg{t});
        calcdev.var_pls2prevpls(t, :) = var(pls2prevpls(t, :));
    end

end

% Display and save data 
pe_asv_4disp = [ave_pls2ent std_pls2ent var_pls2ent];
pi_asv_4disp = [ave_pls2imk std_pls2imk var_pls2imk];
pt_asv_4disp = [ave_pls2tar std_pls2tar var_pls2tar];
tpt_asv_4disp = [ave_dist_pls_tar_off_tar std_dist_pls_tar_off_tar var_dist_pls_tar_off_tar];
pi_pitch_asv_4disp = [ave_pitch std_pitch var_pitch];
pi_roll_asv_4disp = [ave_roll std_roll var_roll];
pi_yaw_asv_4disp = [ave_yaw std_yaw var_yaw];

% Display trigger devations and variance 
disp('===========================================================================');
disp('                       TRIGGER DEVIATIONS & VARIANCE                       ');
disp('===========================================================================');
disp('DISTANCE BETWEEN ENTRY & INSTRUMENT:'); 
disp(imk2ent);
disp('===========================================================================');
disp('DEVIATION BETWEEN TRIGGERED PULSES & ENTRY:'); 
disp('     ave        std        var'); 
disp(pe_asv_4disp); 
disp('===========================================================================');
disp('DEVIATION BETWEEN TRIGGERED PULSES & INSTRUMENT:'); 
disp('     ave        std        var'); 
disp(pi_asv_4disp); 
disp('===========================================================================');
%disp('DISTANCE BETWEEN INSTRUMENT & PREDETERMINED (INTENDED) TARGET:'); 
%disp(imk2tar);
%disp('===========================================================================');
disp('POSITIONAL DEVIATION BETWEEN TRIGGERED PULSES & INSTRUMENT MARKER:'); 
disp('     ave        std        var'); 
disp(pi_asv_4disp);
disp('===========================================================================');
disp('ROTATIONAL DEVIATION (degrees) BETWEEN TRIGGERED PULSES & INSTRUMENT MARKER:');
disp('     ave        std        var'); 
disp(['   ',num2str(pi_pitch_asv_4disp),'       Pitch: tilting the handle forward/backward']);
disp(['  ',num2str(pi_roll_asv_4disp),'       Roll: rotating around the hotspot']);
disp(['  ',num2str(pi_yaw_asv_4disp),'       Yaw: tilting the wings side-to-side']);
disp('AVERAGE COORDINATES OF ESTIMATED STIMULATED SITE (BASED ON PULSE LOCATIONS):'); 
disp('       X        Y        Z'); 
disp(ave_pls_tar);
disp('===========================================================================');
disp('DEVIATION BETWEEN PREDETERMINED (INTENDED) TARGET AND THE ESTIMATED STIMULATED SITE:');
disp('     ave        std        var'); 
disp(tpt_asv_4disp);

names = strrep(trigger_fields, '_', ' ');
names = strrep(names, 'dayMT100', 'STIM');
names = strrep(names, 'dayMT10', 'SHAM');


% Plotting section for distance between pulses and entry/instrument
figure('Color',[1,1,1],'Name', [folder_name, ' - Distance between Pulses and Entry/Instrument'], 'NumberTitle', 'off');
for t = 1:length(to_graph)
    subplot(length(to_graph), 1, t);
    points = 1:numel(pls2ent(to_graph(t), :));
    set(gca, 'YGrid', 'on', 'GridLineStyle' ,'-');
    h = plot(points, pls2ent(to_graph(t),:), 'b'); hold on;
    h = plot(points, pls2imk(to_graph(t),:), 'r'); hold on;
    yline(N, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
    title(names{to_graph(t)}); hold on;
    lh = legend(sprintf('PULSES vs. ENTRY\n %d pulses >%dmm off\n -------------------------', pls2ent_ctOFF(to_graph(t), 1), N),...
                sprintf('PULSES vs. INSTRUMENT\n %d pulses >%dmm off\n -------------------------', pls2imk_ctOFF(to_graph(t), 1), N),...
                sprintf('Threshold Line at %dmm', N));
    set(lh, 'Location', 'eastoutside', 'Orientation', 'vertical');
    ylabel('distance (mm)');
    xlabel(sprintf('pulse count: %s', num2str(numel(trg{to_graph(t)}(:, 1)))));
    xlim([0, numel(trg{to_graph(t)}(:, 1)) + 1]); 
end

% Saving the first figure
savefig(fullfile(sprintf('%s/QAfig_EI', dir_QA)));
saveas(gcf, fullfile(dir_QA, sprintf('QAfig_EI_%s.png', folder_name)));


% Plotting section for distance between target and actual stimulated site
figure('Color',[1,1,1],'Name', [folder_name, ' - Distance between Target and Actual Stimulated Site'], 'NumberTitle', 'off');
for t = 1:length(to_graph)
    subplot(length(to_graph), 1, t);
    points = 1:numel(dist_pls_tar_off_tar(to_graph(t),:));
    set(gca, 'YGrid', 'on', 'GridLineStyle' ,'-');
    h = plot(points, dist_pls_tar_off_tar(to_graph(t),:), 'g'); hold on;
    yline(0, 'm--', 'LineWidth', 1.5);
    yline(T, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
    title(names{to_graph(t)}); hold on;
    lh = legend(sprintf('STIMULATED SITE vs. TARGET\n %d pulses >%dmm off\n -------------------------', dist_pls_tar_off_tar_ctOFF(to_graph(t), 1), T),...
                'Target',...
                sprintf('Threshold Line at %dmm', T));
    set(lh, 'Location', 'eastoutside', 'Orientation', 'vertical');
    ylabel('distance (mm)'); 
    xlabel(sprintf('pulse count: %s', num2str(numel(trg{to_graph(t)}(:, 1)))));
    xlim([0, numel(trg{to_graph(t)}(:, 1)) + 1]); 
end

% Saving the second figure
savefig(fullfile(sprintf('%s/QAfig_STIMSITE', dir_QA)));
saveas(gcf, fullfile(dir_QA, sprintf('QAfig_STIMSITE_%s.png', folder_name)));

% Plotting rotational deviations for each pulse relative to instrument
figure('Color',[1,1,1],'Name', [folder_name, ' - Rotational Deviation between Pulses and Entry/Instrument'], 'NumberTitle', 'off');
for t = 1:length(to_graph)
    subplot(length(to_graph), 1, t);
    points = 1:numel(procrustes_trg2ins{t});
    set(gca, 'YGrid', 'on', 'GridLineStyle' ,'-');
   % h = plot(points, procrustes_trg2ins{to_graph(t)}(:), 'black'); hold on;
    h = plot(points, pitch{to_graph(t)}(:),'m'); hold on;
    h = plot(points, roll{to_graph(t)}(:),'g'); hold on;
    h = plot(points, yaw{to_graph(t)}(:),'b'); hold on;
    yline(0, 'k--', 'LineWidth', 1.5);
   
    title(names{to_graph(t)}); hold on;
    lh = legend({'Pitch (tilt forward/backward)','Roll (rotate around hotspot)','Yaw (tilt side-to-side)'});
    set(lh, 'Location', 'eastoutside', 'Orientation', 'vertical');
    ylabel('deviation (degrees)'); 
    xlabel(sprintf('pulse count: %s', num2str(numel(trg{to_graph(t)}(:, 1)))));
    xlim([0, numel(trg{to_graph(t)}(:, 1)) + 1]); 
end

% Saving the second figure
savefig(fullfile(sprintf('%s/QAfig_ROTATION', dir_QA)));
saveas(gcf, fullfile(dir_QA, sprintf('QAfig_ROTATION_%s.png', folder_name)));



