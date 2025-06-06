%% Load and package passive ephys data

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
    bhv_days = {train_rec_passive.day};
    ephys_days =  bhv_days([train_rec_passive.ephys]);

    data_animal = table;

    for use_rec=1:length(ephys_days)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
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

        % Trial stim values
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));

        % Quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');
 
        % Stim-aligned PSTH
        [~,binned_spikes_stim_align] = ap.psth(spike_times_timelite, ...
            stimOn_times (quiescent_trials), ...
            depth_group);

        % Stim-aligned PSTH (MSN-only)
        msn_spikes = ismember(spike_templates, find(striatum_celltypes.msn));
        [~,binned_msn_spikes_stim_align] = ap.psth(spike_times_timelite(msn_spikes), ...
            stimOn_times(quiescent_trials),  ...
            depth_group(msn_spikes));

        % Stim-aligned average PSTH (each unit)
        stim_align_times = arrayfun(@(x) stimOn_times(trial_stim_values == x & quiescent_trials), ...
            unique(trial_stim_values),'uni',false);
        [unit_event_psths,~] = ap.psth(spike_times_timelite,stim_align_times,spike_templates);

        % Stim-responsive units
        unit_resp_stim_t = [0.05,0.15];

        [~,baseline_spikes] = ap.psth(spike_times_timelite,stim_align_times,spike_templates, ...
            'window',-mean(unit_resp_stim_t),'bin_size',diff(unit_resp_stim_t));
        [~,event_spikes] = ap.psth(spike_times_timelite,stim_align_times,spike_templates, ...
            'window',mean(unit_resp_stim_t),'bin_size',diff(unit_resp_stim_t));

        event_spikes_meandiff = cellfun(@(event,baseline) ...
            squeeze(mean(diff([baseline,event],[],2),1)), ...
            event_spikes,baseline_spikes,'uni',false);

        num_shuffles = 10000;
        event_spikes_meandiff_shuff = cell(size(event_spikes_meandiff,1),num_shuffles);
        for shuffle=1:num_shuffles
            event_spikes_meandiff_shuff(:,shuffle) = cellfun(@(event,baseline) ...
                squeeze(mean(diff(ap.shake([baseline,event],2),[],2),1)), ...
                event_spikes,baseline_spikes,'uni',false);
        end

        unit_resp_p_value = cell(size(event_spikes_meandiff));
        for stim_idx=1:size(event_spikes_meandiff,1)
            curr_rank = tiedrank([event_spikes_meandiff{stim_idx}, ...
                cell2mat(event_spikes_meandiff_shuff(stim_idx,:))]');
            
            unit_resp_p_value{stim_idx} = curr_rank(1,:)'/(num_shuffles+1);
        end

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.trial_stim_values(use_rec) = {trial_stim_values(quiescent_trials)};
        data_animal.depth_group_edges(use_rec) = {depth_group_edges};
        data_animal.unit_depth_group(use_rec) = {unit_depth_group};

        data_animal.binned_spikes_event_align(use_rec) = {binned_spikes_stim_align};
        data_animal.binned_msn_spikes_event_align(use_rec) = {binned_msn_spikes_stim_align};

        data_animal.unit_resp_stim_t(use_rec) = {unit_resp_stim_t};
        data_animal.unit_resp_p_value(use_rec) = {unit_resp_p_value};

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
save_name = fullfile(save_path, 'ephys_passive');
save(save_name, "ephys", "-v7.3");

