%% Load and package task widefield data

%% Set save path and animals

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};


%% Loop through animals, grab and store data

% Initialize cell for data
data_all = cell(length(animals),1);

for animal_idx=1:length(animals)

    animal = animals{animal_idx};

    % Find all task recordings
    disp(['Start ' animal])
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);

    data_animal = table;

    for use_rec=1:length(recordings_task)

        rec_day = recordings_task(use_rec).day;
        rec_time = recordings_task(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording

        % Event-align widefield        
        wf_stim_time = -0.5:0.025:1;

        V_stim_align = interp1(wf_t,wf_V',stimOn_times(1:n_trials) + wf_stim_time);
        V_move_align = interp1(wf_t,wf_V',stim_move_time(1:n_trials) + wf_stim_time);
        V_outcome_align = interp1(wf_t,wf_V',stimOff_times(1:n_trials) + wf_stim_time);

        V_event_align = cat(4,V_stim_align,V_move_align,V_outcome_align);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.trial_outcome(use_rec) = {trial_outcome(1:n_trials)};

        data_animal.wf_stim_time(use_rec) = {wf_stim_time};
        data_animal.V_event_align(use_rec) = {V_event_align};

        disp(['Done day ' num2str(use_rec)])
        
    end
    
    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;
    disp(['Done ' animal])

end

% Concatenate data into one table and save
wf = vertcat(data_all{:});
save_name = fullfile(save_path, 'wf_task');
save(save_name, "wf", "-v7.3");

