%% Load and package behavior data

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

data_all = cell(length(animals), 1);

warning off;
for animal_idx=1:length(animals)
    
    animal = animals{animal_idx};
    disp(animal);

    % Find all task recordings
    workflow = {'stim_wheel_right*'};
    recordings = plab.find_recordings(animal, [], workflow);
    
    data_animal = table;
    for use_rec=1:length(recordings)

        rec_day = recordings(use_rec).day;
        rec_time = recordings(use_rec).recording{end};
        verbose = false;
        load_parts.behavior = true;
        ap.load_recording
    
        % Get association p-value in a few ways
        % (first just grab mean null reaction times for each trial)
        [~,~,~,stim_to_move_nullmean] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'mean');

        % (mean to firstmove)
        [stimwheel_pval_firstmove_mean,stimwheel_rxn_firstmove_mean,stimwheel_rxn_null_firstmove_mean] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'mean');

        % (median to firstmove)
        [stimwheel_pval_firstmove_median,stimwheel_rxn_firstmove_median,stimwheel_rxn_null_firstmove_median] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'median');

        % (mad to firstmove)
        [stimwheel_pval_firstmove_mad,stimwheel_rxn_firstmove_mad,stimwheel_rxn_null_firstmove_mad] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'mad');

        % (mean to lastmove)
        [stimwheel_pval_lastmove_mean,stimwheel_rxn_lastmove_mean,stimwheel_rxn_null_lastmove_mean] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_lastmove, 'mean');

        % (mad to lastmove)
        [stimwheel_pval_lastmove_mad,stimwheel_rxn_lastmove_mad,stimwheel_rxn_null_lastmove_mad] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_lastmove, 'mad');

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        % (firstmove mean stats)
        data_animal.stimwheel_pval_firstmove_mean(use_rec) = {stimwheel_pval_firstmove_mean};
        data_animal.stimwheel_rxn_firstmove_mean(use_rec) = {stimwheel_rxn_firstmove_mean};
        data_animal.stimwheel_rxn_null_firstmove_mean(use_rec) = {stimwheel_rxn_null_firstmove_mean};

        % (firstmove median stats)
        data_animal.stimwheel_pval_firstmove_median(use_rec) = {stimwheel_pval_firstmove_median};
        data_animal.stimwheel_rxn_firstmove_median(use_rec) = {stimwheel_rxn_firstmove_median};
        data_animal.stimwheel_rxn_null_firstmove_median(use_rec) = {stimwheel_rxn_null_firstmove_median};

        % (firstmove mad stats)
        data_animal.stimwheel_pval_firstmove_mad(use_rec) = {stimwheel_pval_firstmove_mad};
        data_animal.stimwheel_rxn_firstmove_mad(use_rec) = {stimwheel_rxn_firstmove_mad};
        data_animal.stimwheel_rxn_null_firstmove_mad(use_rec) = {stimwheel_rxn_null_firstmove_mad};

        % (lastmove mean stats)
        data_animal.stimwheel_pval_lastmove_mean(use_rec) = {stimwheel_pval_lastmove_mean};
        data_animal.stimwheel_rxn_lastmove_mean(use_rec) = {stimwheel_rxn_lastmove_mean};
        data_animal.stimwheel_rxn_null_lastmove_mean(use_rec) = {stimwheel_rxn_null_lastmove_mean};

        % (lastmove mad stats)
        data_animal.stimwheel_pval_lastmove_mad(use_rec) = {stimwheel_pval_lastmove_mad};
        data_animal.stimwheel_rxn_lastmove_mad(use_rec) = {stimwheel_rxn_lastmove_mad};
        data_animal.stimwheel_rxn_null_lastmove_mad(use_rec) = {stimwheel_rxn_null_lastmove_mad};

        % (reaction times)
        data_animal.stim_to_move(use_rec) = {stim_to_move(1:n_trials)};
        data_animal.stim_to_outcome(use_rec) = {stim_to_outcome(1:n_trials)};
        data_animal.trial_outcome(use_rec) = {logical(trial_outcome(1:n_trials))};

        % (null reaction times)
        data_animal.stim_to_move_nullmean(use_rec) = {stim_to_move_nullmean(1:n_trials)};
        
        % Print progress
        ap.print_progress_fraction(use_rec,length(recordings))
    end

    data_all{animal_idx} = data_animal;

end
warning on;

bhv = vertcat(data_all{:});

save_name = fullfile(save_path, 'bhv');
save(save_name, "bhv", "-v7.3");

fprintf('Saved %s\n',save_name);

