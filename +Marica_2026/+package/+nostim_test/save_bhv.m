%% Load and package behavior data

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data\nostim';

animals = { ...
    'AP024','AP026'};

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
    
        % Get association p-value 
        % (for no stim: this is reaction time with vs without stim)
        trial_opacity = logical(vertcat(trial_events.values(1:n_trials).TrialOpacity));

        rxn_stat = diff(ap.groupfun(@mean,stim_to_move(1:n_trials),trial_opacity));
        n_shuff = 10000;
        rxn_null_stat_distribution = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            rxn_null_stat_distribution(curr_shuff) = ...
                diff(ap.groupfun(@mean,stim_to_move,ap.shake(trial_opacity)));
        end

        rxn_stat_rank = tiedrank(vertcat(rxn_stat,rxn_null_stat_distribution));
        rxn_stat_p = rxn_stat_rank(1)./(n_shuff+1);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        % (association stat)
        data_animal.trial_opacity(use_rec) = {trial_opacity};
        data_animal.stimwheel_pval(use_rec) = rxn_stat_p;
      
        % (reaction times)
        data_animal.stim_to_move(use_rec) = {stim_to_move(1:n_trials)};
        data_animal.stim_to_outcome(use_rec) = {stim_to_outcome(1:n_trials)};
        data_animal.trial_outcome(use_rec) = {logical(trial_outcome(1:n_trials))};

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

