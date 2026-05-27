%% Load and package ephys properties

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
        load_parts.behavior = true;
        load_parts.ephys = true;
        ap.load_recording

        % Get striatum boundaries
        AP_longstriatum_find_striatum_depth
        mua_length = 200;
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(end);

        if any(isnan(striatum_depth))
            % Skip recording if no striatum detected
            warning(['No striatum: ' animal ' ' rec_day])
            continue
        end

        depth_group = discretize(spike_depths,depth_group_edges);
        unit_depth_group = discretize(template_depths, depth_group_edges);

        % Single unit labels
        single_unit_idx = strcmp(template_qc_labels, 'singleunit');

        % Classify striatal celltype
        AP_longstriatum_classify_striatal_units

        % Grab pre-saved ACGs
        [~, ~, acg] = bc.ep.loadSavedProperties(qMetrics_path);
        acg_matrix = table2array(acg(good_templates,:));

        % Store cell type properties
        data_animal.waveform(use_rec) = {waveforms};
        data_animal.acg(use_rec) = {acg_matrix};
        data_animal.waveformDuration_peakTrough_us(use_rec) = {ephysProperties.waveformDuration_peakTrough_us};
        data_animal.postSpikeSuppression_ms(use_rec) = {ephysProperties.postSpikeSuppression_ms};
        data_animal.mean_firingRate(use_rec) = {ephysProperties.mean_firingRate};

        data_animal.striatal_units(use_rec) = {striatal_units};
        data_animal.single_unit_idx(use_rec) = {single_unit_idx};
        data_animal.str_tan_idx(use_rec) = {striatum_celltypes.tan};
        data_animal.str_fsi_idx(use_rec) = {striatum_celltypes.fsi};
        data_animal.str_msn_idx(use_rec) = {striatum_celltypes.msn};

        fprintf('Done: %s day %d/%d\n',animal,use_rec,length(recordings_task));
    end

    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;

end

% Concatenate data into one table and save
ephys_properties = vertcat(data_all{:});
save_name = fullfile(save_path, 'ephys_properties');
save(save_name, "ephys_properties", "-v7.3");
fprintf('Saved: %s\n',save_name);
