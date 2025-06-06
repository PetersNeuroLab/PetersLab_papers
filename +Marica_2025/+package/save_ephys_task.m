%% Load and package task ephys data

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

data_all = cell(length(animals),1);

for animal_idx=1:length(animals)
    
    animal = animals{animal_idx};

    % Find all task recordings
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);

    data_animal = table;

    for use_rec=1:length(recordings_task)

        rec_day = recordings_task(use_rec).day;
        rec_time = recordings_task(use_rec).recording{end};
        verbose = true;
        load_parts.behavior = true;
        load_parts.ephys = true;
        ap.load_recording

        % Get striatum boundaries
        AP_longstriatum_find_striatum_depth
        mua_length = 200;
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(end);
        
        if any(isnan(striatum_depth))
            % Skip recording if no striatum detected
            warning(['Undefined str depth ' animal ' ' rec_day])
            data_animal.animal(use_rec) = {animal};
            data_animal.rec_day(use_rec) = {rec_day};
            continue
        end

        depth_group = discretize(spike_depths,depth_group_edges);
        unit_depth_group = discretize(template_depths, depth_group_edges);

        % Single unit labels
        single_unit_idx = strcmp(template_qc_labels, 'singleunit'); 

        % Classify striatal celltype
        AP_longstriatum_classify_striatal_units

        % Event-aligned PSTHs
        [~,binned_spikes_stim_align] = ap.psth(spike_times_timelite, ...
            stimOn_times(1:n_trials),  ...
            depth_group);

        [~,binned_spikes_move_align] = ap.psth(spike_times_timelite, ...
            stim_move_time(1:n_trials),  ...
            depth_group);

        [~,binned_spikes_outcome_align] = ap.psth(spike_times_timelite, ...
            stimOff_times(1:n_trials),  ...
            depth_group);

        binned_spikes_event_align = cat(4,binned_spikes_stim_align,binned_spikes_move_align,binned_spikes_outcome_align);

        % Event-aligned PSTHs (MSN-only)
        msn_spikes = ismember(spike_templates, find(striatum_celltypes.msn));

         [~,binned_msn_spikes_stim_align] = ap.psth(spike_times_timelite(msn_spikes), ...
            stimOn_times(1:n_trials),  ...
            depth_group(msn_spikes));

        [~,binned_msn_spikes_move_align] = ap.psth(spike_times_timelite(msn_spikes), ...
            stim_move_time(1:n_trials),  ...
            depth_group(msn_spikes));

        [~,binned_msn_spikes_outcome_align] = ap.psth(spike_times_timelite(msn_spikes), ...
            stimOff_times(1:n_trials),  ...
            depth_group(msn_spikes));

        binned_msn_spikes_event_align = cat(4,binned_msn_spikes_stim_align,binned_msn_spikes_move_align,binned_msn_spikes_outcome_align);

        % Event-aligned average PSTHs (each unit)
        [unit_stim_psths,~] = ap.psth(spike_times_timelite,stimOn_times(1:n_trials),spike_templates);
        [unit_move_psths,~] = ap.psth(spike_times_timelite,stim_move_time(1:n_trials),spike_templates);
        [unit_outcome_psths,~] = ap.psth(spike_times_timelite,stimOff_times(1:n_trials),spike_templates);

        unit_event_psths = cat(3,unit_stim_psths,unit_move_psths,unit_outcome_psths);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.depth_group_edges(use_rec) = {depth_group_edges};
        data_animal.unit_depth_group(use_rec) = {unit_depth_group};

        data_animal.binned_spikes_event_align(use_rec) = {binned_spikes_event_align};
        data_animal.binned_msn_spikes_event_align(use_rec) = {binned_msn_spikes_event_align};
        data_animal.unit_event_psths(use_rec) = {unit_event_psths};

        data_animal.single_unit_idx(use_rec) = {single_unit_idx}; 
        data_animal.str_tan_idx(use_rec) = {striatum_celltypes.tan}; 
        data_animal.str_fsi_idx(use_rec) = {striatum_celltypes.fsi}; 
        data_animal.str_msn_idx(use_rec) = {striatum_celltypes.msn}; 

        disp(['Done day ' num2str(use_rec)])

    end

    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;
    disp(['Done ' animal]);

end

% Concatenate data into one table and save
ephys = vertcat(data_all{:});
save_name = fullfile(save_path, 'ephys_task');
save(save_name, "ephys", "-v7.3");

