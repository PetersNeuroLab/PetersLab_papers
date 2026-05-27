%% Load and package passive widefield data

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

data_all = cell(length(animals),1);

for animal_idx=1:length(animals)

    animal = animals{animal_idx};

    % Find passive recording days that also have task
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);

    data_animal = table;

    for use_rec=1:length(train_rec_passive)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        load_parts.behavior = true;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording

        % Trial stim values
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));
        
        % Quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % Event-align widefield        
        wf_stim_time = -0.5:0.025:1;
        V_stim_align = interp1(wf_t,wf_V',stimOn_times(quiescent_trials) + wf_stim_time);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.trial_stim_values(use_rec) = {trial_stim_values(quiescent_trials)};

        data_animal.wf_stim_time(use_rec) = {wf_stim_time};
        data_animal.V_event_align(use_rec) = {V_stim_align};

        disp(['Done day ' num2str(use_rec)])
        
    end
    
    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;
    disp(['Done ' animal])

end

% Concatenate data into one table and save
wf = vertcat(data_all{:});
save_name = fullfile(save_path, 'wf_passive');
save(save_name, "wf", "-v7.3");

