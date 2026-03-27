%% Load and package cortex-striatum maps

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
        ap.load_recording;
  
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

        % Get time bins corresponding to widefield frame exposures
        % (skip the beginning and end of the recording to avoid artifacts)
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        % Bin spikes in depth and time
        binned_spikes = zeros(max(depth_group),length(time_bins)-1);
        for curr_depth = 1:max(depth_group)
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        %% Regress ctx fluorescence to striatal MUA

        % Normalize MUA by standard devation
        % (to have kernels in comparable units)
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

        % Set parameters for regression
        use_svs = 1:200; 
        kernel_t = [-0.2,0.2]; 
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 10;
        zs = [false,false];
        cvfold = 1;
        use_constant = true;
        return_constant = false;

        % Interpolate the widefield to the center of each time bin to match the
        % MUA timepoints
        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Do regression from cortex to striatum
        [cortex_kernel,predicted_spikes,explained_var] = ...
            ap.regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % Convert kernel V into pixels
        cortex_kernel_px = squeeze(plab.wf.svd2px(wf_U(:,:,use_svs),cortex_kernel));

        % Define map as max kernel across lags
        cortex_striatum_map = permute(max(cortex_kernel_px,[],3),[1,2,4,3]);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.depth_group_edges(use_rec) = {depth_group_edges};
        data_animal.cortex_striatum_map(use_rec) = {cortex_striatum_map}; 
        data_animal.explained_var(use_rec) = {explained_var.total};

        disp(['Done day ' num2str(use_rec)])
        
    end

    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;
    disp(['Done ' animal])

end

% Concatenate data into one table and save
ctx_str_maps = vertcat(data_all{:});
save_name = fullfile(save_path, 'ctx_str_maps');
save(save_name, "ctx_str_maps", "-v7.3");

