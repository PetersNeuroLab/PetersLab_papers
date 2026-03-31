%% Generate figures for Marica et al 2026
%
% Can save figures and stats if flags are checked.
%
% Data packaged in +Marica_2026.package namespace
%
% Loads/preps relevant data at the start of each figure.
% Doesn't re-load data if already loaded.
% 
% Each cell can be run independently, or whole script can be run through.
% If script run through, figure/stat saving is set at top.
%
% Dependencies are in: 
% petersaj/AP_scripts_peterslab
% PetersNeuroLab/PetersLab_analysis

tic

%% Set path to save figures and print stats (if script thru-run)

% Path to save figures
fig_savepath = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2026','figures','matlab_figs');

% Filename to print stats
stat_savefn = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2026','figures','stats.txt');

% Set flag to overwrite save (restrict to AP username)
if strcmp(getenv('USERNAME'),'petersa')

    % (safety catch: save only on AP's computer)
    fig_overwrite_confirm = strcmp(questdlg('Overwite saved figures?', ...
        'Confirm save','No','Yes','No'),'Yes');
    if fig_overwrite_confirm
        % (check or make directory) 
        if ~exist(fig_savepath,'dir')
            mkdir(fig_savepath);
        end
        % (turn on flag to save figs)
        fig_save_flag = true;
        % (set function to save figures)
        save_figs = @() arrayfun(@(curr_fig) saveas(curr_fig, ...
            fullfile(fig_savepath,strrep(curr_fig.Name,' ','_')),'fig'), ...
            findall(0,'Type','figure'));
    end

    stat_overwrite_confirm = strcmp(questdlg('Overwrite stats?', ...
        'Confirm save','No','Yes','No'),'Yes');
    if stat_overwrite_confirm
        % (create stat file for writing)
        stat_fid = fopen(stat_savefn,'w');
        % (set stats to print to stat file)
        print_stat = @(varargin) fprintf(stat_fid,varargin{:});
    end

end


%% [Fig 1B-C; Fig S1B] Reaction mean/index, learning day histogram

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Plot reaction time and association index, split within day
n_daysplit = 3;
use_rxn = cellfun(@(x) ~isnan(x),bhv.stim_to_move_nullmean,'uni',false);

rxn_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,x(idx), ...
    ap.quantile_bin(sum(idx),n_daysplit)),bhv.stim_to_move,use_rxn,'uni',false)')';

rxn_null_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,x(idx), ...
    ap.quantile_bin(sum(idx),n_daysplit)),bhv.stim_to_move_nullmean,use_rxn,'uni',false)')';

rxn_idx_daysplit = (rxn_null_mean_daysplit-rxn_mean_daysplit)./ ...
    (rxn_null_mean_daysplit+rxn_mean_daysplit);

[rxn_daysplit_mean,rxn_group_x] = ap.groupfun(@mean,rxn_mean_daysplit, ...
    bhv.days_from_learning);
rxn_null_daysplit_mean = ap.groupfun(@mean,rxn_null_mean_daysplit, ...
    bhv.days_from_learning);

rxn_idx_daysplit_mean = ap.groupfun(@mean,rxn_idx_daysplit, ...
    bhv.days_from_learning);
rxn_idx_daysplit_sem = ap.groupfun(@AP_sem,rxn_idx_daysplit, ...
    bhv.days_from_learning);

plot_days = -3:2;
plot_day_idx = ismember(rxn_group_x,plot_days);

figure('Name','Fig 1 rxn'); tiledlayout(2,1);
rxn_group_x_daysplit = rxn_group_x+(0:n_daysplit)./n_daysplit;

nexttile; hold on; set(gca,'YScale','log');
plot(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
plot(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_null_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1),'r','linewidth',2);
xline(0,'r');
ylabel('Reaction time');
xlabel('Day from learning');

nexttile;
errorbar(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_idx_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1), ...
    reshape(padarray(rxn_idx_daysplit_sem(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
xline(0,'r');
ylabel('Association index');
xlabel('Day from learning');
ap.prettyfig;


% Plot histogram of learning days
n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
figure('Name','Fig S1 hist');
histogram(n_learned_day,-0.5:max(n_learned_day)+0.5,'FaceColor','k','EdgeColor','none');
ylabel('Number of mice');
xlabel('Days to learn');
ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 1--\n');

stat_day = -1;
stat_use_data = bhv.days_from_learning == stat_day;
p = anovan(reshape(rxn_idx_daysplit(stat_use_data,:),[],1), ...
    reshape(repmat(1:3,sum(stat_use_data),1),[],1),'display','off'); 
print_stat('Rxn idx 1-way ANOVA day %d, p = %.2g\n',stat_day,p);

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 1D] Example widefield/units

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Choose animal and day to plot
use_animal = 'AP023';

animal_days = find(strcmp(bhv.animal,use_animal));
plot_days = [2,length(animal_days)];

figure('Name','Fig 1 examples');
h = tiledlayout(2,3);
title(h,sprintf('%s, days %d,%d',use_animal,plot_days));
for curr_rec_idx = plot_days

    % Load data from example recording
    use_rec = animal_days(curr_rec_idx);
    animal = bhv.animal{use_rec};
    rec_day = bhv.rec_day{use_rec};
    recordings = plab.find_recordings(animal,rec_day,'stim_wheel*');
    rec_time = recordings.recording{end};
    load_parts.ephys = true;
    ap.load_recording;

    % Plot widefield average
    nexttile(h);

    wf_day_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
    mean_image_fn = fullfile(wf_day_path,sprintf('meanImage_blue.npy'));
    wf_avg = plab.wf.wf_align(readNPY(mean_image_fn),animal,rec_day);
    imagesc(wf_avg);
    ap.wf_draw('cortex','y');
    axis image off;
    clim([0,12000]);

    % Plot average domain map colored and combined
    domain_avg = ap.groupfun(@mean,ctx_str_maps.cortex_striatum_map{use_rec},[],[],domain_idx_rec{use_rec});

    domain_color = {'180','60','310'};
    col_lim = [0,0.01];
    domain_colored = nan([size(domain_avg),3]);
    for curr_domain = 1:n_domains
        curr_colormap = ap.colormap(['W',domain_color{curr_domain}],[],2);
        curr_map_gray = 1+round(mat2gray(domain_avg(:,:,curr_domain),col_lim).*(size(curr_colormap,1)-1));
        domain_colored(:,:,curr_domain,:) = reshape(curr_colormap(curr_map_gray,:),size(curr_map_gray,1),size(curr_map_gray,2),3,[]);
    end

    domain_colored_combined = squeeze(min(domain_colored,[],3));
    
    nexttile(h);
    image(domain_colored_combined)
    axis image off;
    ap.wf_draw('cortex',[0.5,0.5,0.5]);

    % Plot units and overlay clustering
    domain_color_rgb = [1,0,0;0,1,0;0,0,1];
    
    ax = nexttile(h); hold on;

    curr_domain_idx = domain_idx_rec{use_rec};
    domain_im = permute(domain_color_rgb(curr_domain_idx,:),[1,3,2]);
    imagesc(ax,[],ctx_str_maps.depth_group_edges{use_rec}/1000,domain_im);
    ax.YDir = 'reverse';

    spike_rate = (accumarray(findgroups(spike_templates),1)+1)/ ...
        range(spike_times_timelite);
    unit_dots = scatter( ...
        spike_rate,template_depths(unique(spike_templates))/1000,20,'k','filled');
    ylabel('Depth (mm)')
    xlabel('Spike rate')
    set(gca,'XScale','log');

end

ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 1E] Corticostriatal map pre/post learning

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Make map indicies (not done in load data)
wf_map_grp = struct;
wf_map_grp.ld = cell2mat(cellfun(@(idx,maps) repmat(idx,size(maps,3).*~isempty(maps),1), ...
    num2cell(bhv.days_from_learning),ctx_str_maps.cortex_striatum_map,'uni',false));
wf_map_grp.animal = cell2mat(cellfun(@(idx,maps) repmat(idx,size(maps,3).*~isempty(maps),1), ...
    num2cell(grp2idx(bhv.animal)),ctx_str_maps.cortex_striatum_map,'uni',false));

% Group pre/post learn days
plot_day_bins = [-Inf,-0,Inf];
plot_day_grp = discretize(max(wf_map_grp.ld,-inf),plot_day_bins);

% Get and plot average maps
[wf_map_avg,wf_map_avg_grp] = ap.nestgroupfun({@mean,@mean}, ...
    reshape(cat(3,ctx_str_maps.cortex_striatum_map{:}),prod(U_size),[])', ...
    wf_map_grp.animal,[plot_day_grp,domain_idx]);

figure('Name','Fig 1 cstr maps');
domain_color = {'180','60','310'};
h = tiledlayout(n_domains,length(plot_day_bins)-1,'TileSpacing','none');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        curr_data_idx = ismember(wf_map_avg_grp,[curr_day,curr_domain],'rows');
        
        nexttile;
        imagesc(reshape(wf_map_avg(curr_data_idx,:),size(U_master,[1,2])));
        axis image off;
        colormap(gca,ap.colormap(['W',domain_color{curr_domain}],[],2));
        clim([0,0.01]);
        ap.wf_draw('cortex',[0.5,0.5,0.5]);
    end
end

ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 1--\n');

% Average maps in domain within recording
[maps_animal,maps_grp] = ap.groupfun(@mean, ...
    cat(3,ctx_str_maps.cortex_striatum_map{:}), ...
    [],[],[wf_map_grp.animal,wf_map_grp.ld,domain_idx,plot_day_grp]);

% Median filter and rectify > 0.5x max weight
maps_animal_filt = medfilt3(maps_animal,[3,3,1]);
maps_animal_filt(maps_animal_filt < max(maps_animal_filt,[],[1,2]).*0.5) = 0;

% Get grid of correlations
maps_corr_grid = corrcoef(reshape(maps_animal_filt,[],size(maps_animal_filt,3)));
animal_grid = repmat(maps_grp(:,1),1,size(maps_grp,1));

% Loop through domains, get corr across day vs day-shuffle w/i animal
print_stat('Map correlation across day group vs shuffle:\n');
for curr_domain = 1:n_domains

    curr_idx = ...
        maps_grp(:,1) == maps_grp(:,1)' & ...   % same animal
        (maps_grp(:,3) == curr_domain & ...
        maps_grp(:,3)' == curr_domain) & ...    % in domain
        maps_grp(:,4) ~= maps_grp(:,4)' & ...   % across day groups
        tril(true(size(maps_grp,1)),-1);        % single comparisons

    map_corr = mean(ap.groupfun(@mean,maps_corr_grid(curr_idx),animal_grid(curr_idx)));

    [~,~,shuff_grp] = unique(maps_grp(:,[1,3]),'rows');
    n_shuff = 10000;
    map_corr_shuff = nan(n_shuff,1);
    for curr_shuff = 1:n_shuff
        day_grp_shuff = ap.shake(maps_grp(:,4),1,shuff_grp);

        curr_idx = ...
            maps_grp(:,1) == maps_grp(:,1)' & ...   % same animal
            (maps_grp(:,3) == curr_domain & ...
            maps_grp(:,3)' == curr_domain) & ...    % in domain
            day_grp_shuff ~= day_grp_shuff' & ...   % across day groups
            tril(true(size(maps_grp,1)),-1);        % single comparisons

        map_corr_shuff(curr_shuff) = mean(ap.groupfun(@mean, ...
            maps_corr_grid(curr_idx),animal_grid(curr_idx)));
    end

    stat_rank = tiedrank([map_corr;map_corr_shuff]);
    stat_p = stat_rank(1)/(n_shuff+1);
    print_stat('Domain %d, p = %.2g\n',curr_domain,stat_p);

end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 1F] Striatal MUA pre/post learning

%%% Load data for figure
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
learn_colormap = ap.colormap('BKR',3);
prepost_colormap = max(0,learn_colormap([1,end],:)-0.2);

% Plot average task activity (pre/post learning)
align_t = { ...
    [-0.2,0.15]; % stim
    [-0.05,0.4]; % move
    [-0.1,0.5]}; % reward

figure('Name','Fig 1 mua'); h = tiledlayout(n_domains,3,'TileSpacing','tight');
for curr_depth = 1:n_domains

    curr_trials = striatum_mua_grp.domain_idx == curr_depth;

    for curr_align = 1:3
        use_t = isbetween(psth_t,align_t{curr_align}(1),align_t{curr_align}(2));
        nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
        curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
            striatum_mua(curr_trials,use_t,curr_align), ...
            striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
        curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
            striatum_mua(curr_trials,use_t,curr_align), ...
            striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
        ap.errorfill(psth_t(use_t),curr_data_mean',curr_data_sem');
        xline(0);
    end

end
% (link all y, and x of same-alignment)
linkaxes(h.Children,'y');
for ax = 1:3
    linkaxes(h.Children(ax:3:end),'x');
end
% (set same data aspect to have same x-size)
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));
ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 1--\n');
print_stat('MUA 2-way ANOVA, group x time interaction:\n')
for curr_align = 1:3

    use_t = isbetween(psth_t,align_t{curr_align}(1),align_t{curr_align}(2));

    [curr_data_animalmean,curr_data_animalmean_grp] = ...
        ap.groupfun(@nanmean,striatum_mua(:,use_t,curr_align), ...
        [striatum_mua_grp.animal,plot_day_grp,striatum_mua_grp.domain_idx]);

    for curr_domain = 1:n_domains
        curr_grps_use = curr_data_animalmean_grp(:,3) == curr_domain;
        stat_grp_day = repmat(curr_data_animalmean_grp(curr_grps_use,2),1,sum(use_t));
        stat_grp_t = repmat(1:sum(use_t),size(stat_grp_day,1),1);

        p = anovan(reshape(curr_data_animalmean(curr_grps_use,:),[],1), ...
            {stat_grp_t(:),stat_grp_day(:)},'model','interaction','display','off');
        
        print_stat('Align %d, domain %d, p = %.2g\n',curr_align,curr_domain,p(3));
    end

end


% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end

%% [Diagram] Cortex ROIs

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

figure('Name','Diagram cortex rois');tiledlayout(n_domains,1,'tilespacing','none')
for curr_domain = 1:n_domains
    nexttile;
    imagesc(striatum_wf_roi(:,:,curr_domain));
    axis image off; ap.wf_draw('ccf',[0.5,0.5,0.5]);
    colormap(ap.colormap('WG'));
end
ap.prettyfig;


% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 2A,D; Fig S6A] Task cortex + striatum trial heatmap

%%% Load data for figure
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);


% Plot heatmaps sorted by reaction times
heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)

figure('Name','Fig 2 heatmaps'); 
h = tiledlayout(n_domains*2,max(cortex_plot_day_grp), ...
    'TileIndexing','ColumnMajor','TileSpacing','Tight');
for curr_day = unique(cortex_plot_day_grp)'
    for curr_domain = 1:n_domains

        % Cortex
        nexttile;
        curr_trials = find(cortex_plot_day_grp == curr_day);
       
        [sorted_rxn,sort_idx] = sort(wf_grp.rxn(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_domain,1),heatmap_smooth));
        colormap(gca,ap.colormap('WG'));
        clim([0,1e-2]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end

        % Striatum
        nexttile;
        curr_trials = find(striatum_plot_day_grp == curr_day & striatum_mua_grp.domain_idx == curr_domain);
      
        [sorted_rxn,sort_idx] = sort(striatum_mua_grp.rxn(curr_trials));      
        imagesc(psth_t,[],movmean(striatum_mua(curr_trials(sort_idx),:,1),heatmap_smooth));
        colormap(gca,ap.colormap('WK'));
        clim([0,3]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');

    end
end
linkaxes(h.Children,'x');
xlim(h.Children,[-0.5,1])
ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 2B,E; Fig S6B] Task cortex + striatum average PSTH w/o movement

%%% Load data for figure
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% NaN-out activity after movement onset (minus leeway time)
move_leeway = 0.1; % time pre-movement to exclude
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

% Get average and SEM no-movement activity
[striatum_mua_nomove_avg,striatum_mua_nomove_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
striatum_mua_nomove_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_nomove_avg,wf_striatum_roi_nomove_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);
wf_striatum_roi_nomove_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);

str_day_color = ap.colormap('KW',length(plot_day_bins)-1);
ctx_day_color = ap.colormap('KG',length(plot_day_bins)-1);
figure('Name','Fig 2 psth'); h = tiledlayout(n_domains*2,1,'TileSpacing','tight');
for curr_domain = 1:n_domains

    nexttile; hold on; set(gca,'ColorOrder',ctx_day_color);
    ap.errorfill(wf_t,wf_striatum_roi_nomove_avg(:,:,curr_domain)', ...
        wf_striatum_roi_nomove_sem(:,:,curr_domain)');

    nexttile; hold on; set(gca,'ColorOrder',str_day_color);
    plot_data_idx = striatum_mua_nomove_avg_grp(:,2) == curr_domain;
    ap.errorfill(psth_t,striatum_mua_nomove_avg(plot_data_idx,:)', ...
        striatum_mua_nomove_sem(plot_data_idx,:)');

end

linkaxes(h.Children(1:2:end));
linkaxes(h.Children(2:2:end));
xlim(h.Children,[-0.1,0.5]);
ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 2--\n');
for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

    stim_t = [0,0.2];
    cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
    striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

    [wf_striatum_roi_nomove_animal,wf_striatum_roi_nomove_animal_grp] = ...
        ap.groupfun(@nanmean,wf_striatum_roi_nomove, ...
        [wf_grp.animal,cortex_plot_day_grp]);
    wf_striatum_roi_nomove_animal_tmax = max(wf_striatum_roi_nomove_animal(:,cortex_stim_t,:),[],2);

    cortex_stat_usedata = ismember(wf_striatum_roi_nomove_animal_grp(:,2),compare_day_grps);
    cortex_stat_meas = permute(diff(ap.groupfun(@mean,wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:), ...
        wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,2)),[],1),[3,2,1]);

    [striatum_mua_nomove_animal,striatum_mua_nomove_animal_grp] = ...
        ap.groupfun(@nanmean,striatum_mua_nomove, ...
        [striatum_mua_grp.animal,striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
    [~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[1,3]),'rows');
    striatum_mua_nomove_animal_tmax = max(striatum_mua_nomove_animal(:,striatum_stim_t),[],2);

    striatum_stat_usedata = ismember(striatum_mua_nomove_animal_grp(:,2),compare_day_grps);
    striatum_stat_meas = ap.nestgroupfun({@mean,@diff},striatum_mua_nomove_animal_tmax(striatum_stat_usedata), ...
        striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

    n_shuff = 10000;
    cortex_stat_null = nan(n_domains,n_shuff);
    striatum_stat_null = nan(n_domains,n_shuff);
    for curr_shuff = 1:n_shuff

        curr_ctx_shuff = ...
            ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1, ...
            wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,1));
        cortex_stat_null(:,curr_shuff) = ...
            permute(diff(ap.groupfun(@mean,curr_ctx_shuff, ...
            wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,2)),[],1),[3,2,1]);     

        curr_str_shuff = ...
            ap.shake(striatum_mua_nomove_animal_tmax(striatum_stat_usedata),1, ...
            striatum_shuffgroup(striatum_stat_usedata));
        striatum_stat_null(:,curr_shuff) = ...
            ap.nestgroupfun({@mean,@diff},curr_str_shuff, ...
            striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

    end

    print_stat('Average PSTH tmax, day grps %d,%d:\n',compare_day_grps);

    cortex_stat_rank = tiedrank([cortex_stat_meas,cortex_stat_null]')';
    cortex_stat_p = 1-cortex_stat_rank(:,1)/(n_shuff+1);
    cortex_stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],{'','*'});
    for curr_domain = 1:n_domains
        print_stat('CTX%d p = %.2g%s\n',curr_domain,cortex_stat_p(curr_domain),cortex_stat_sig{curr_domain});
    end

    striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
    striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);
    striatum_stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],{'','*'});
    for curr_domain = 1:n_domains
        print_stat('STR%d p = %.2g%s\n',curr_domain,striatum_stat_p(curr_domain),striatum_stat_sig{curr_domain});
    end

end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 2C,F; Fig S6C] Task cortex + striatum max, split within session

%%% Load data for figure
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-2:0,0.5];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 4;

striatum_split_idx = cell2mat(cellfun(@(n_trials,n_mua) ...
    reshape(repmat(ap.quantile_bin(n_trials,n_split),1,n_mua),[],1), ...
    num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

cortex_split_idx = cell2mat(cellfun(@(n_trials) ...
    ap.quantile_bin(n_trials,n_split), ...
    num2cell(cortex_trials_rec_n),'uni',false));

% For bins that contain multiple days: don't split within day
% (doesn't makes sense to average splits across days)
striatum_split_idx(ismember(striatum_plot_day_grp,find(diff(plot_day_bins) > 1))) = 1;
cortex_split_idx(ismember(cortex_plot_day_grp,find(diff(plot_day_bins) > 1))) = 1;

% Get activity: trial mean > time max > animal mean 
stim_t = [0,0.2];

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > rxn_cutoff;
[striatum_activity_mean,striatum_activity_mean_grp] = ap.groupfun(@mean,striatum_mua(striatum_use_trials,:,1), ...
    [striatum_mua_grp.animal(striatum_use_trials),striatum_plot_day_grp(striatum_use_trials), ...
    striatum_split_idx(striatum_use_trials),striatum_mua_grp.domain_idx(striatum_use_trials)]);

striatum_activity_max = max(striatum_activity_mean(:,isbetween(psth_t,stim_t(1),stim_t(2))),[],2);
[striatum_activity_max_avg,striatum_activity_max_grp] = ap.nestgroupfun({@mean,@mean},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));
striatum_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));

% (cortex)
cortex_use_trials = wf_grp.rxn > rxn_cutoff;
[cortex_activity_mean,cortex_activity_mean_grp] = ap.groupfun(@mean,wf_striatum_roi(cortex_use_trials,:,:,1), ...
    [wf_grp.animal(cortex_use_trials),cortex_plot_day_grp(cortex_use_trials),cortex_split_idx(cortex_use_trials)]);

cortex_activity_max = permute(max(cortex_activity_mean(:,isbetween(wf_t,stim_t(1),stim_t(2)),:),[],2),[1,3,2]);

[cortex_activity_max_avg,cortex_activity_max_grp] = ap.nestgroupfun({@mean,@mean},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));
cortex_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));

% Plot cortex and striatum overlaid
figure('Name','Fig 2 session split');
h = tiledlayout(n_domains,length(plot_day_bins)-1,'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;

        yyaxis left;
        curr_cortex_data_idx = cortex_activity_max_grp(:,1) == curr_day;
        if sum(curr_cortex_data_idx) > 1
            ap.errorfill(1:sum(curr_cortex_data_idx),cortex_activity_max_avg(curr_cortex_data_idx,curr_domain), ...
                cortex_activity_max_sem(curr_cortex_data_idx),[0,0.6,0],0.5,false,2);
        else
            errorbar(cortex_activity_max_avg(curr_cortex_data_idx,curr_domain), ...
                cortex_activity_max_sem(curr_cortex_data_idx),'color',[0,0.6,0],'linewidth',4, ...
                'CapSize',0);
        end

        yyaxis right;
        curr_striatum_data_idx = striatum_activity_max_grp(:,1) == curr_day & ...
            striatum_activity_max_grp(:,3) == curr_domain;
        errorbar(striatum_activity_max_avg(curr_striatum_data_idx), ...
            striatum_activity_max_sem(curr_striatum_data_idx),'k','linewidth',2, ...
            'CapSize',0);

        if curr_domain == 1
            title(sprintf('%d:%g',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end

% Set axis limits
for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'left');
end
linkaxes(h.Children);
ylim(h(1).Children,[0,8e-3]);
[h.Children.YColor] = deal([0,0.5,0]);

for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'right');
end
linkaxes(h.Children);
ylim(h(1).Children,[-0.2,8]);
[h.Children.YColor] = deal([0,0,0]);

xlim(h(1).Children,[0.8,n_split+0.2])

ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 2--\n');
print_stat('Session-split 1-way ANOVA:\n');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        stat_data_idx = ismember(striatum_activity_mean_grp(:,[2,4]),[curr_day,curr_domain],'rows');
        p = anovan(striatum_activity_max(stat_data_idx),striatum_activity_mean_grp(stat_data_idx,3),'display','off');
        stat_sig = discretize(p < 0.05,[0,1,Inf],["","*"]);
        print_stat('STR %d, day grp %d, p = %.2g%s\n',curr_domain,curr_day,p,stat_sig);
    end
end
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        stat_data_idx = cortex_activity_mean_grp(:,2) == curr_day;
        p = anovan(cortex_activity_max(stat_data_idx,curr_domain),cortex_activity_mean_grp(stat_data_idx,3),'display','off');
        stat_sig = discretize(p < 0.05,[0,1,Inf],["","*"]);
        print_stat('CTX %d, day grp %d, p = %.2g%s\n',curr_domain,curr_day,p,stat_sig);
    end
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 3B,C,E,F; Fig S6D; Fig S7] Passive Cortex + striatum PSTH

%%% Load data for figure
load_dataset = 'passive';
Marica_2026.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
striatum_plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

[striatum_mua_avg,striatum_mua_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);
striatum_mua_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);
wf_striatum_roi_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

unique_stim = unique(striatum_mua_grp.stim);
stim_color = {'KB','KW','KR'};

figure('Name','Fig 3 psth');
h = tiledlayout(n_domains*2,max(striatum_plot_days_grp)*length(unique_stim), ...
    'TileIndexing','ColumnMajor','TileSpacing','tight');
for curr_stim = unique_stim'

    [~,curr_stim_idx] = ismember(curr_stim,unique_stim);
    day_colormap = ap.colormap(stim_color{curr_stim_idx},length(plot_day_bins)-1);

    for curr_day = 1:length(plot_day_bins)-1
        for curr_domain = 1:n_domains

            % (cortex)
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_day & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;
            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_domain), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_domain), ...
                day_colormap(curr_day,:));

            % (striatum)
            nexttile; axis off;
            curr_data_idx = striatum_mua_avg_grp(:,1) == curr_day & ...
                striatum_mua_avg_grp(:,2) == curr_stim & ...
                striatum_mua_avg_grp(:,3) == curr_domain;
            ap.errorfill(psth_t,striatum_mua_avg(curr_data_idx,:),striatum_mua_sem(curr_data_idx,:), ...
                day_colormap(curr_day,:));

        end
    end
end
linkaxes(h.Children(1:2:end),'y');
linkaxes(h.Children(2:2:end),'y');
xlim(h.Children,[0,0.5]);
ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 3G-H; Fig S6E] Passive cortex + striatum max

%%% Load data for figure
load_dataset = 'passive';
Marica_2026.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

% (cortex)
[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    (1:size(wf_striatum_roi,1))', ...
    [wf_grp.animal,cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,cortex_stim_t,:),[],2);

[wf_striatum_roi_max_avg,wf_striatum_roi_max_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

wf_striatum_roi_max_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

% (striatum)
[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@mean,striatum_mua, ...
    [striatum_mua_grp.animal,striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

striatum_mua_dayavg_tmax = max(striatum_mua_dayavg(:,striatum_stim_t),[],2);

[striatum_mua_max_avg,striatum_mua_max_grp] = ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));
striatum_mua_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));

% Plot activity by day
figure('Name','Fig 3 tmax');
h = tiledlayout(n_domains,2,'TileSpacing','tight');
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = 1:size(striatum_wf_roi,3)
    % (cortex)
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    errorbar(binned_days_x, ...
        reshape(wf_striatum_roi_max_avg(:,:,curr_domain),[],length(plot_day_bins)-1)', ...
        reshape(wf_striatum_roi_max_sem(:,:,curr_domain),[],length(plot_day_bins)-1)', ...
        'linewidth',2);
    ylabel('Max \DeltaF/F_0');
    axis padded
    xline(0);

    % (striatum)
    nexttile; hold on; 
    set(gca,'ColorOrder',stim_colormap);
    curr_str_idx = striatum_mua_max_grp(:,3) == curr_domain;
    errorbar(binned_days_x, ...
        reshape(striatum_mua_max_avg(curr_str_idx),[],length(plot_day_bins)-1)', ...
        reshape(striatum_mua_max_sem(curr_str_idx),[],length(plot_day_bins)-1)', ...
        'linewidth',2);
    ylabel('Max \DeltaR/R_0');
    axis padded
    xline(0);
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy');

ap.prettyfig;

% ~~~ STATS ~~~
% (compare day i to i+1)
print_stat('\n--FIG 3--\n');
print_stat('PSTH max\n');
for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

    cortex_stat_usedata = ismember(wf_striatum_roi_grp(:,2),compare_day_grps);
    [cortex_stat_meas,cortex_stat_grp] = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
        wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

    striatum_stat_usedata = ismember(striatum_mua_dayavg_grp(:,2),compare_day_grps);
    [striatum_stat_meas,striatum_stat_grp] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata), ...
        striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:end));

    [~,~,cortex_shuff_grp] = unique(wf_striatum_roi_grp(cortex_stat_usedata,[1,3]),'rows');
    [~,~,striatum_shuff_grp] = unique(striatum_mua_dayavg_grp(striatum_stat_usedata,[1,3,4]),'rows');

    n_shuff = 10000;
    cortex_stat_null = nan(size(cortex_stat_meas,1),n_shuff,n_domains);
    striatum_stat_null = nan(length(striatum_stat_meas),n_shuff);
    for curr_shuff = 1:n_shuff
        cortex_data_shuff = ap.shake(wf_striatum_roi_max(cortex_stat_usedata,:,:),1,cortex_shuff_grp);
        cortex_stat_null(:,curr_shuff,:) = ap.nestgroupfun({@mean,@diff},cortex_data_shuff, ...
            wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

        striatum_data_shuff = ap.shake(striatum_mua_dayavg_tmax(striatum_stat_usedata),1,striatum_shuff_grp);
        striatum_stat_null(:,curr_shuff) = ap.nestgroupfun({@mean,@diff},striatum_data_shuff, ...
            striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:end));
    end

    cortex_stat_rank = permute(tiedrank(permute([cortex_stat_meas,cortex_stat_null],[2,1,3])),[2,1,3]);
    cortex_stat_p = 1-cortex_stat_rank(:,1,:)/(n_shuff+1);

    striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
    striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);

    print_stat('Cortex: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(cortex_stat_grp(:,1))'
            curr_stat_idx = ismember(cortex_stat_grp,curr_stim,'rows');
            print_stat('ROI%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,cortex_stat_grp(curr_stat_idx),cortex_stat_p(curr_stat_idx,1,curr_domain),stat_sig(curr_stat_idx,1,curr_domain));
        end
    end

    print_stat('Striatum: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(striatum_stat_grp(:,1))'
            curr_stat_idx = ismember(striatum_stat_grp,[curr_stim,curr_domain],'rows');
            print_stat('D%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,striatum_stat_grp(curr_stat_idx,1),striatum_stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
        end
    end
end

% (compare C to R stim within-day)
print_stat('PSTH max stim comparison\n');
compare_stim = [0,90];

cortex_stat_usedata = ismember(wf_striatum_roi_grp(:,3),compare_stim);
[cortex_stat_meas,cortex_stat_grp] = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
    wf_striatum_roi_grp(cortex_stat_usedata,3),wf_striatum_roi_grp(cortex_stat_usedata,2));

striatum_stat_usedata = ismember(striatum_mua_dayavg_grp(:,3),compare_stim);
[striatum_stat_meas,striatum_stat_grp] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata), ...
    striatum_mua_dayavg_grp(striatum_stat_usedata,3),striatum_mua_dayavg_grp(striatum_stat_usedata,[2,4]));

[~,~,cortex_shuff_grp] = unique(wf_striatum_roi_grp(cortex_stat_usedata,[1,2]),'rows');
[~,~,striatum_shuff_grp] = unique(striatum_mua_dayavg_grp(striatum_stat_usedata,[1,2,4]),'rows');

n_shuff = 10000;
cortex_stat_null = nan(size(cortex_stat_meas,1),n_shuff,n_domains);
striatum_stat_null = nan(length(striatum_stat_meas),n_shuff);
for curr_shuff = 1:n_shuff
    cortex_data_shuff = ap.shake(wf_striatum_roi_max(cortex_stat_usedata,:,:),1,cortex_shuff_grp);
    cortex_stat_null(:,curr_shuff,:) = ap.nestgroupfun({@mean,@diff},cortex_data_shuff, ...
    wf_striatum_roi_grp(cortex_stat_usedata,3),wf_striatum_roi_grp(cortex_stat_usedata,2));

    striatum_data_shuff = ap.shake(striatum_mua_dayavg_tmax(striatum_stat_usedata),1,striatum_shuff_grp);
    striatum_stat_null(:,curr_shuff) = ap.nestgroupfun({@mean,@diff},striatum_data_shuff, ...
    striatum_mua_dayavg_grp(striatum_stat_usedata,3),striatum_mua_dayavg_grp(striatum_stat_usedata,[2,4]));
end

cortex_stat_rank = permute(tiedrank(permute([cortex_stat_meas,cortex_stat_null],[2,1,3])),[2,1,3]);
cortex_stat_p = 1-cortex_stat_rank(:,1,:)/(n_shuff+1);

striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);

print_stat('Cortex: Stim %d vs %d\n',compare_stim);
stat_sig = discretize(cortex_stat_p,[-Inf,0.05,0.95,Inf],["*-","","*+"]);
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        curr_stat_idx = ismember(cortex_stat_grp,curr_day,'rows');
        print_stat('ROI%d, Day %d, p = %.2g%s\n', ...
            curr_domain,cortex_stat_grp(curr_stat_idx),cortex_stat_p(curr_stat_idx,1,curr_domain),stat_sig(curr_stat_idx,1,curr_domain));
    end
end

print_stat('Striatum: Stim %d vs %d\n',compare_stim);
stat_sig = discretize(striatum_stat_p,[-Inf,0.05,0.95,Inf],["*-","","*+"]);
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        curr_stat_idx = ismember(striatum_stat_grp,[curr_day,curr_domain],'rows');
        print_stat('D%d, Day %d, p = %.2g%s\n', ...
            curr_domain,striatum_stat_grp(curr_stat_idx,1),striatum_stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
    end
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end

%% [Fig 4A] Task/passive pre-learn PSTHs

% Get task traces (NaN-out activity after movement onset)
load_dataset = 'task';
Marica_2026.figures.load_data;

% Set day binning 
plot_day_bins = [-2,-1];

striatum_plot_day_grp = discretize(striatum_mua_grp.ld,plot_day_bins);
cortex_plot_day_grp = discretize(wf_grp.ld,plot_day_bins);

move_leeway = 0.1; % time pre-movement to exclude
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

[striatum_mua_task_avg,striatum_mua_task_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
striatum_mua_task_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_task_avg,wf_striatum_roi_task_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);
wf_striatum_roi_task_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);

% Get passive traces
% (turn on retain data to load task+passive)
load_dataset_retain = true;
load_dataset = 'passive';
Marica_2026.figures.load_data;
load_dataset_retain = false;

striatum_plot_day_grp = discretize(striatum_mua_grp.ld,plot_day_bins);
cortex_plot_day_grp = discretize(wf_grp.ld,plot_day_bins);

[striatum_mua_passive_avg,striatum_mua_passive_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_day_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);
striatum_mua_passive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_day_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_passive_avg,wf_striatum_roi_passive_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);
wf_striatum_roi_passive_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

% Plot task/passive overlay
plot_passive_stim = 90;
figure('Name','Fig 4 task passive psths');
tiledlayout(n_domains,2,'TileSpacing','compact');
for curr_domain = 1:n_domains

    str_task_idx = striatum_mua_task_avg_grp(:,2) == curr_domain;

    wf_passive_idx = wf_striatum_roi_passive_avg_grp(:,2) == plot_passive_stim;
    str_passive_idx = ismember(striatum_mua_passive_avg_grp(:,2:3),[plot_passive_stim,curr_domain],'rows');

    nexttile; hold on;
    ap.errorfill(wf_t,wf_striatum_roi_task_avg(:,:,curr_domain),wf_striatum_roi_task_sem(:,:,curr_domain),'k');
    ap.errorfill(wf_t,wf_striatum_roi_passive_avg(wf_passive_idx,:,curr_domain), ...
        wf_striatum_roi_passive_sem(wf_passive_idx,:,curr_domain),'r');
    title(sprintf('Cortex %d',curr_domain));

    nexttile; hold on;
    ap.errorfill(psth_t,striatum_mua_task_avg(str_task_idx,:),striatum_mua_task_sem(str_task_idx,:),'k');
    ap.errorfill(psth_t,striatum_mua_passive_avg(str_passive_idx,:), ...
        striatum_mua_passive_sem(str_passive_idx,:),'r');
    title(sprintf('Striatum %d',curr_domain));
end

legend({'Task','','Passive'});

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end

%% [Fig 4B-D - Load/prep data for below cells]

% Grab day-binned task/passive stim responses

% ~~ Set up data structure params
load_dataset = 'task';
Marica_2026.figures.load_data;

% Set up parameters for activity grids [animal x day x domain x stim]
data_grid_params = struct;

data_grid_params.stim_t = [0,0.2];
data_grid_params.cortex_stim_t = isbetween(wf_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.striatum_stim_t = isbetween(psth_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.rxn_cutoff = 0.3;

% data_grid_params.ld_bins = [-Inf,-2:0,1,Inf];
data_grid_params.ld_bins = [-Inf,-2:0,0]; % end at LD0 - task after has few trials
% data_grid_params.ld_bins = [-Inf,-2,0,2,Inf]; % same as paper - best
% % data_grid_params.ld_bins = [-Inf,-2:2,2];
data_grid_params.ld_unique = 1:(length(data_grid_params.ld_bins)-1);
data_grid_params.grid_size = [length(unique(bhv.animal)),length(data_grid_params.ld_unique),n_domains];

% Set up data structure
data_grids = struct;

% ~~ Get task activity

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > data_grid_params.rxn_cutoff;
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials & ~isnan(striatum_ld_idx))));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_task = accumarray(striatum_rec_grp,striatum_rec_tmax,data_grid_params.grid_size,[],NaN);

% (widefield)
wf_use_trials = wf_grp.rxn > data_grid_params.rxn_cutoff;
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx].* ...
    ap.nanout(~(wf_use_trials & ~isnan(wf_ld_idx))));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_task = cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),data_grid_params.grid_size(1:2),[],NaN('single')),1:n_domains,'uni',false),[1,3,2]));

% ~~ Get passive activity
% (manually clear workspace to keep previously loaded)
clearvars -except  ...
    load_dataset fig_save_flag stat_fid print_stat save_figs ... % (standard in load_data)
    data_grid_params data_grids                                  % (task data from above)

% (turn on retain data to load task+passive)
load_dataset_retain = true;
load_dataset = 'passive';
Marica_2026.figures.load_data;
load_dataset_retain = false;

% (striatum)
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);
[~,striatum_stim_idx] = ismember(striatum_mua_grp.stim,unique(striatum_mua_grp.stim));

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx,striatum_stim_idx].* ...
    ap.nanout(isnan(striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_passive = accumarray(striatum_rec_grp,striatum_rec_tmax,[data_grid_params.grid_size,max(striatum_stim_idx)],[],NaN);

% (widefield)
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);
[~,wf_stim_idx] = ismember(wf_grp.stim,unique(wf_grp.stim));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx,wf_stim_idx].* ...
    ap.nanout(isnan(wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_passive = permute(cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),[data_grid_params.grid_size(1:2),length(unique(wf_grp.stim))], ...
    [],NaN('single')),1:n_domains,'uni',false),[1,3,4,2])),[1,2,4,3]);

% ~~ Get widefield stim kernels
U_master = plab.wf.load_master_U;
load(fullfile(data_path,'wf_kernels'));

n_vs = size(wf_kernels.task_kernels{1},1);

wf_grid_idx = [grp2idx(bhv.animal),discretize(bhv.days_from_learning,data_grid_params.ld_bins)];
wf_grid_idx_use = ~any(isnan(wf_grid_idx),2);

% (task)
wf_kernel_roi_task = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[3,2,1]), ...
    wf_kernels.task_kernels,'uni',false));

wf_kernel_roi_task_tmax = permute(max(wf_kernel_roi_task,[],2),[1,3,2]);

data_grids.wf_kernel_roi_task = ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_task_tmax(wf_grid_idx_use,domain),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2]));

% (passive)
wf_kernel_roi_passive = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[4,2,1,3]), ...
    wf_kernels.passive_kernels,'uni',false));
wf_kernel_roi_passive_tmax = permute(max(wf_kernel_roi_passive,[],2),[1,3,4,2]);

data_grids.wf_kernel_roi_passive = cell2mat(permute(arrayfun(@(stim) ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_passive_tmax(wf_grid_idx_use,domain,stim),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2])),1:size(wf_kernel_roi_passive_tmax,3),'uni',false), [1,3,4,2]));


%% |--> [Fig 4B] Activity by context

plot_str = [1,2];
plot_ctx = 2;

area_labels= ["Striatum " + string(num2cell(plot_str)), ...
    "Cortex " + string(num2cell(plot_ctx))];

use_stim = 3;
plot_passive_data = cat(3,data_grids.striatum_passive(:,:,plot_str,use_stim), ...
    data_grids.wf_roi_passive(:,:,plot_ctx,use_stim));
plot_task_data = cat(3,data_grids.striatum_task(:,:,plot_str), ...
    data_grids.wf_roi_task(:,:,plot_ctx));

figure('Name','Fig 4 context activity'); hold on; 
xlim([0,1]);ylim(xlim);
n_areas = length([plot_str,plot_ctx]);
area_colors = ap.colormap('tube',n_areas);
task_passive_scale = cell(n_areas,1);
h_dot = gobjects(n_areas,1);
for curr_area = 1:n_areas
        % Max-normalize
        curr_act_norm = rescale(cat(3, ...
            plot_passive_data(:,:,curr_area),plot_task_data(:,:,curr_area)));
        
        % Fit scale passive to task
        %  (all data)
        nonan_idx = ~any(isnan(curr_act_norm),3);
        curr_act_norm_vec = reshape(curr_act_norm(repmat(nonan_idx,1,1,2)),[],2);
        curr_scale = curr_act_norm_vec(:,1)\curr_act_norm_vec(:,2);

        % (by animal)
        curr_scale_animal = nan(size(curr_act_norm,1),1);
        for curr_animal = 1:size(curr_act_norm,1)
            nonan_idx = ~any(isnan(curr_act_norm(curr_animal,:,:)),3);
            if sum(nonan_idx) < 2
                continue
            end
            curr_scale_animal(curr_animal) = ...
                curr_act_norm(curr_animal,nonan_idx,1)'\ ...
                curr_act_norm(curr_animal,nonan_idx,2)';
        end
        task_passive_scale{curr_area} = curr_scale_animal;
        
        % Scatter and plot slope fit with error shading
        h_dot(curr_area) = scatter(reshape(curr_act_norm(:,:,1),[],1), ...
            reshape(curr_act_norm(:,:,2),[],1),30, ...
            area_colors(curr_area,:),'filled');

        slope_mean = nanmean(task_passive_scale{curr_area});
        slope_sem = ap.sem(task_passive_scale{curr_area});

        % (plot average slope across animals)
        % ap.errorfill([0,1],[0,slope_mean],[0,slope_sem],area_colors(curr_area,:));

        % (plot slope fit across all data)
        h_line = refline(curr_scale);
        h_line.Color = area_colors(curr_area,:);
        h_line.LineWidth = 2;

        line([0,1],[0,1],'color',[0.5,0.5,0.5]);
        xlim(ylim);
        axis square; xlim([0,1]);ylim([0,1]);
        xlabel('Passive');ylabel('Task');
end
legend(h_dot,area_labels,'location','se');
ap.prettyfig;

% Plot slopes by region
figure('Name','Fig 4 context slope'); hold on;

swarmchart(cell2mat(cellfun(@(x,y) repelem(x,length(y),1), ...
    num2cell(1:length(task_passive_scale))',task_passive_scale,'uni',false)), ...
    cell2mat(task_passive_scale),'filled','MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(cellfun(@nanmean,task_passive_scale(:)), ...
    cellfun(@ap.sem,task_passive_scale(:)),'k','linewidth',2,'CapSize',0);

set(gca,'XTick',1:n_areas,'XTickLabels',area_labels);

yline(1);
ylabel('Task/passive scale');
ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG 4--\n');

sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

for curr_domain = 1:n_areas-1
    stat_p = signrank(task_passive_scale{end},task_passive_scale{curr_domain});
    print_stat('Signrank context mPFC vs Striatum %d: p = %.2g%s\n', ...
        curr_domain,stat_p,sig_flag(stat_p));
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% |--> [Fig 4C-D] mPFC vs striatum by context

% Normalize task/passive separately
% (normalize data to task LD 0)
norm_bin = find(data_grid_params.ld_bins>=0,1);

str_task_norm = data_grids.striatum_task./data_grids.striatum_task(:,norm_bin,:);
ctx_task_norm = data_grids.wf_roi_task./data_grids.wf_roi_task(:,norm_bin,:);

use_stim = 3;
str_passive_norm = data_grids.striatum_passive(:,:,:,use_stim)./data_grids.striatum_task(:,norm_bin,:);
ctx_passive_norm = data_grids.wf_roi_passive(:,:,:,use_stim)./data_grids.wf_roi_task(:,norm_bin,:);

str_passive_norm = data_grids.striatum_passive(:,:,:,use_stim)./(data_grids.striatum_passive(:,norm_bin,:,3)+1);
ctx_passive_norm = data_grids.wf_roi_passive(:,:,:,use_stim)./(data_grids.wf_roi_passive(:,norm_bin,:,3)+0.001);

% Plot mPFC vs striatum 1/2
figure('Name','Fig 4 mpfc v striatum'); tiledlayout(1,2);
plot_wf_v_str = @(wf_data,str_data,col) ...
    errorbar(nanmean(str_data,1),nanmean(wf_data,1), ...
    ap.sem(wf_data,1),ap.sem(wf_data,1),ap.sem(str_data,1),ap.sem(str_data,1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30,'capsize',0);
plot_wf = 2;
str_col = {'k','r'};
for curr_str = 1:2
    nexttile(1); hold on;
    plot_wf_v_str(ctx_task_norm(:,:,plot_wf),str_task_norm(:,:,curr_str),str_col{curr_str});
    title('Task');
    nexttile(2); hold on;
    plot_wf_v_str(ctx_passive_norm(:,:,plot_wf),str_passive_norm(:,:,curr_str),str_col{curr_str});
    title('Passive');
end
for curr_tile = 1:2
    nexttile(curr_tile); axis square; line(xlim,xlim,'color',[0.5,0.5,0.5])
    xlabel('Striatum');ylabel('Cortex');
    legend("Striatum "+string(num2cell(1:2)));
end
ap.prettyfig;

% ~~~ STATS ~~~
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

% Compare mPFC to striatum 1/2
print_stat('\n--FIG 4--\n');

print_stat('2-way ANOVA cortex v striatum:\n');

use_ctx = 2;
use_str = [1,2];

[~,group_ld] = ndgrid(1:data_grid_params.grid_size(1), ...
    1:data_grid_params.grid_size(2));

for curr_context = ["Task","Passive"]
    for curr_str = use_str
        switch curr_context
            case "Task"
                curr_ctx_data = ctx_task_norm(:,:,use_ctx);
                curr_str_data = str_task_norm(:,:,use_str(curr_str));
            case "Passive"
                curr_ctx_data = ctx_passive_norm(:,:,use_ctx);
                curr_str_data = str_passive_norm(:,:,use_str(curr_str));
        end

        use_data = ~any(isnan(cat(3,curr_ctx_data,curr_str_data)),3);

        p = anovan([curr_ctx_data(use_data);curr_str_data(use_data)], ...
            [reshape(ones(sum(use_data,'all'),2).*[1:2],[],1), ...
            repmat(group_ld(use_data),2,1)], ...
            'display','off');

        print_stat('%s ctx %d v. str %d: group p = %.2g%s, day p = %.2g%s\n', ...
            curr_context,use_ctx,curr_str,p(1),sig_flag(p(1)),p(2),sig_flag(p(2)))
    end
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig 5] Passive striatum cell types and single-unit responses

%%% Load data for figure
load_dataset = 'passive';
Marica_2026.figures.load_data;
%%%

plot_domains = 1:2;

% Set stim to compare, and colors
compare_stim = [2,3]; % C, R stim
stim_1_color = [0,0,0.8];
stim_2_color = [0.8,0,0];

% Get max response in window
use_t = psth_t > 0.05 & psth_t < 0.15;
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.99;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Set celltypes to loop through
striatum_celltypes = ["msn","fsi","tan"];

% Scatter plots R vs C (days combined)
figure('Name','Fig 5 scatter');
h = tiledlayout(length(striatum_celltypes),3, ...
    'TileSpacing','tight');

example_units = nan(length(striatum_celltypes),2);
h_scatter = gobjects(length(striatum_celltypes),3);
for curr_celltype = striatum_celltypes

    curr_units = ...
        find(ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        any(striatum_units_responsive(:,compare_stim),2) & ...
        striatum_sua_grp.(curr_celltype));

    curr_dot_size = 15;
    curr_dot_color = ...
        stim_1_color.*striatum_units_responsive(curr_units,compare_stim(1)) + ...
        stim_2_color.*striatum_units_responsive(curr_units,compare_stim(2));

    [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);
    
    % (split plots by color for grouping in fig prep)
    [~,plot_ax] = ismember(curr_dot_color, ...
        [stim_1_color;stim_2_color;stim_1_color+stim_2_color],'rows');

    for curr_ax = 1:3
        h_scatter(curr_celltype_idx,curr_ax) = nexttile; hold on;
        plot_units = curr_units(plot_ax == curr_ax);

        scatter(unit_psth_max(plot_units,2),unit_psth_max(plot_units,3), ...
            curr_dot_size,curr_dot_color(plot_ax == curr_ax,:),'filled');

        axis square;
        if curr_ax > 1
            h_scatter(curr_celltype_idx,curr_ax).XTick = '';
            h_scatter(curr_celltype_idx,curr_ax).YTick = '';
        else
            title(curr_celltype);
        end
    end
end
for curr_celltype = 1:length(striatum_celltypes)
    linkaxes(h_scatter(curr_celltype,:));
    ax_range = [0,max(cell2mat(ylim(h_scatter(curr_celltype,:))),[],'all')];
    [h_scatter(curr_celltype,:).XLim] = deal(ax_range);
    [h_scatter(curr_celltype,:).YLim] = deal(ax_range);
end
ap.prettyfig;


% Fraction [R,R+C] as stacked barplot over days
figure('Name','Fig 5 bar');
h = tiledlayout(length(striatum_celltypes),1, ...
    'TileSpacing','tight');

for curr_celltype = striatum_celltypes

    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    striatum_units_responsive_2_12 = permute(all(striatum_units_responsive(:,compare_stim) == ...
        cat(3,[0,1],[1,1]),2),[1,3,2]);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive_2_12(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    % Error bars: combine R and C+R-responsive
    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +any(striatum_units_responsive_2_12(use_units,:),2), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    nexttile; hold on;
    set(gca,'ColorOrder',[stim_2_color;stim_1_color+stim_2_color]);
    bar(unit_responsive_mean,'stacked');
    errorbar(sum(unit_responsive_mean,2),unit_responsive_sem, ...
        'k','linestyle','none','linewidth',2);
    ylabel('Frac. responsive units');
    legend(["R","C+R"]);
    title(curr_celltype);

end
ap.prettyfig;

% Plot example PSTHs

% (grab example units to plot)
plot_unit_prctile = 80;
for curr_celltype = striatum_celltypes
    for curr_stim = compare_stim
        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            striatum_units_responsive(:,curr_stim) & ...
            ~striatum_units_responsive(:,setxor(curr_stim,compare_stim)) & ...
            striatum_sua_grp.(curr_celltype);

        curr_unit_idx = find(curr_units);
        [~,sort_idx] = sort(unit_psth_max(curr_units,curr_stim));

        curr_example_unit = curr_unit_idx(sort_idx(ceil(sum(curr_units).*plot_unit_prctile/100)));

        % (circle example units on scatter)
        [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);

        scatter(h_scatter(curr_celltype_idx,1), ...
            unit_psth_max(curr_example_unit,2),unit_psth_max(curr_example_unit,3), ...
            45,[0.5,0.5,0.5],'linewidth',2);

        [~,curr_stim_idx] = ismember(curr_stim,compare_stim);
        example_units(curr_celltype_idx,curr_stim_idx) = curr_example_unit;
    end
end

% (load and plot example units)
figure('Name',sprintf('Fig 5 %dprctile',plot_unit_prctile));
h_units = tiledlayout(1,numel(example_units));
for curr_unit = reshape(example_units',1,[])
    animal = ephys.animal{striatum_sua_grp.rec(curr_unit)};
    rec_day = ephys.rec_day{striatum_sua_grp.rec(curr_unit)};
    unit_id = striatum_sua_grp.unit_id(curr_unit);

    task_flag = true; % (plot task activity)
    quiescence_flag = false; % (plot all trials in passive);
    Marica_2026.figures.psth_figure(animal,rec_day,unit_id,task_flag,quiescence_flag,h_units);    
end
xlim(vertcat(h_units.Children.Children),[-0.1,0.5])
unit_plots = {h_units.Children.Children};
raster_plots = cellfun(@(x) x(1:end-1),unit_plots,'uni',false);
linkaxes(vertcat(raster_plots{:}),'y');
ap.prettyfig;


% ~~~ STATS ~~~
n_shuff = 10000;
print_stat('\n--FIG 5--\n');
for curr_celltype = striatum_celltypes

    % Compare R+C overlap to shuffling R/C responsiveness
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    [~,~,overlap_shuff_grp] = unique([ ...
        striatum_sua_grp.animal, ...
        unit_plot_day_grp, ...
        striatum_sua_grp.domain_idx, ...
        any(striatum_units_responsive(:,compare_stim),2) & ...
        striatum_sua_grp.(curr_celltype)],'rows');

    responsive_overlap_meas = mean(ap.groupfun(@mean, ...
        +(all(striatum_units_responsive(use_units,compare_stim),2)), ...
        striatum_sua_grp.animal(use_units)));

    responsive_overlap_shuff = nan(n_shuff,1);
    for curr_shuff = 1:n_shuff
        striatum_units_responsive_shuff = ap.shake(striatum_units_responsive(use_units,:),1,overlap_shuff_grp(use_units));
        responsive_overlap_shuff(curr_shuff) = mean(ap.groupfun(@mean, ...
            +(all(striatum_units_responsive_shuff(:,compare_stim),2)), ...
            striatum_sua_grp.animal(use_units)));
    end

    stat_rank = tiedrank([responsive_overlap_meas;responsive_overlap_shuff]);
    stat_p = stat_rank(1)/(n_shuff+1);

    stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    print_stat('%s R+C overlap vs shuffle: p = %.2g%s\n',curr_celltype, ...
        stat_p,stat_sig);

    % Compare R fraction (R-only and R+C) across days
    for curr_compare_day = 1:length(plot_day_bins)-2

        compare_day_grps = curr_compare_day+[0,1];

        use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            ismember(unit_plot_day_grp,compare_day_grps) & ...
            striatum_sua_grp.(curr_celltype);

        [~,~,r_shuff_grp] = unique([striatum_sua_grp.animal, ...
            striatum_sua_grp.domain_idx, ...
            striatum_sua_grp.(curr_celltype)],'rows');

        r_frac_meas = diff(ap.nestgroupfun({@mean,@mean}, ...
            +(all(striatum_units_responsive(use_units,3),2)), ...
            striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units)));

        r_frac_shuff = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            unit_plot_day_grp_shuff = ap.shake(unit_plot_day_grp(use_units,:),1,r_shuff_grp(use_units));
            r_frac_shuff(curr_shuff) = diff(ap.nestgroupfun({@mean,@mean}, ...
                +(all(striatum_units_responsive(use_units,3),2)), ...
                striatum_sua_grp.animal(use_units),unit_plot_day_grp_shuff));
        end

        stat_rank = tiedrank([r_frac_meas;r_frac_shuff]);
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
        print_stat('%s R-frac day %d vs %d: p = %.2g%s\n',curr_celltype,compare_day_grps,stat_p,stat_sig);
    end

end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S1A] Example wheel velocity

animal = 'AM021';

use_workflow = 'stim_wheel*';
recordings = plab.find_recordings(animal,[],use_workflow);

plot_days = [2,length(recordings)];

surround_time = [-2,2];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

figure('Name','Fig S1 example wheel');
h = tiledlayout(2,2,'TileSpacing','tight');
for curr_day = plot_days

    % Load data
    rec_day = recordings(curr_day).day;
    rec_time = recordings(curr_day).recording{end};
    load_parts = struct;
    load_parts.behavior = true;
    ap.load_recording;

    % Plot example wheel velocity
    plot_t = [220,285]; % s, example time to plot
    plot_t_idx = isbetween(timelite.timestamps,plot_t(1),plot_t(2));

    nexttile;
    plot(timelite.timestamps(plot_t_idx),wheel_velocity(plot_t_idx),'k','linewidth',1);
    xline(stimOn_times(isbetween(stimOn_times,plot_t(1),plot_t(2))),'g');
    xline(stimOff_times(isbetween(stimOff_times,plot_t(1),plot_t(2))),'r');
    xline(reward_times(isbetween(reward_times,plot_t(1),plot_t(2))),'b');
        
    % Plot stim-aligned wheel velocity
    align_times = stimOn_times;
    pull_times = align_times + surround_time_points;
    event_aligned_wheel_vel = interp1(timelite.timestamps,wheel_velocity,pull_times);

    nexttile;
    imagesc(surround_time_points,[],event_aligned_wheel_vel);
    xline(0,'k');
    clim(4500.*[-1,1])
    colormap(gca,ap.colormap('BWR'));

end

linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy');

ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S1C] Correct performance rate by day

%%% Load non-activity data
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Plot reaction time and association index, split within day
n_daysplit = 3;

outcome_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,+x, ...
    ap.quantile_bin(length(x),n_daysplit)),bhv.trial_outcome,'uni',false)')';

[outcome_daysplit_mean,outcome_group_x] = ...
    ap.groupfun(@mean,outcome_mean_daysplit,bhv.days_from_learning);
outcome_daysplit_sem = ap.groupfun(@AP_sem,outcome_mean_daysplit, ...
    bhv.days_from_learning);

plot_days = -3:2;
plot_day_idx = ismember(outcome_group_x,plot_days);

figure('Name','Fig S1 performance');
outcome_group_x_daysplit = outcome_group_x+(0:n_daysplit)./n_daysplit;
errorbar(reshape(outcome_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(outcome_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1), ...
    reshape(padarray(outcome_daysplit_sem(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
xline(0,'r');
ylabel('Fraction correct');
xlabel('Day from learning');
ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S1D-E] Rate of non-stim movements and P(stim|move)

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

use_stat = 'mean';

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = 'stim_wheel_right_stage\d';
    recordings = plab.find_recordings(animal,[],use_workflow);

    n_bins = 3;

    move_rate = nan(length(recordings),n_bins);
    move_rate_stim = nan(length(recordings),n_bins);
    move_rate_nonstim = nan(length(recordings),n_bins);
    p_cue_given_move = nan(length(recordings),n_bins);
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get P(cue|wheel_move)
        move_onsets = timelite.timestamps(diff(wheel_move) == 1);
        move_onset_prev_stim = interp1(stimOn_times,stimOn_times,move_onsets,'previous','extrap');
        move_prevstim_t = move_onsets - move_onset_prev_stim;

        move_poststim = move_prevstim_t < 0.4;

        % Bin by time within session
        session_bin_edges = linspace(timelite.timestamps(1),timelite.timestamps(end),n_bins+1);
        move_session_bins = discretize(move_onsets,session_bin_edges);

        % Get probability of stim around movements
        p_cue_given_move(curr_recording,:) = ap.groupfun(@mean,+move_poststim,move_session_bins);

        % Get rate of movement with/without stimulus present       
        [move_onset_tidx_valid,move_onset_tidx] = ismember(move_onsets,timelite.timestamps);
        if any(~move_onset_tidx_valid)
            error('Non-valid move onset time index');
        end

        stim_move_idx = photodiode_bw_interp(move_onset_tidx)==1;      
        t_session_bins = discretize(timelite.timestamps,session_bin_edges);

        % % (rate relative to stim on/off period)
        % move_rate_stim(curr_recording,:) = ap.groupfun(@sum,+stim_move_idx,move_session_bins)./ ...
        %     (ap.groupfun(@nanmean,+(photodiode_bw_interp==1),t_session_bins).*diff(session_bin_edges)');
        % move_rate_nonstim(curr_recording,:) = ap.groupfun(@sum,+~stim_move_idx,move_session_bins)./ ...
        %     (ap.groupfun(@nanmean,+(photodiode_bw_interp==0),t_session_bins).*diff(session_bin_edges)');

        % (rate relative to entire time bin e.g. like absolute number)
        move_rate(curr_recording,:) = ap.groupfun(@length,move_onsets,move_session_bins)./diff(session_bin_edges)';
        move_rate_stim(curr_recording,:) = ap.groupfun(@sum,+stim_move_idx,move_session_bins)./diff(session_bin_edges)';
        move_rate_nonstim(curr_recording,:) = ap.groupfun(@sum,+~stim_move_idx,move_session_bins)./diff(session_bin_edges)';

        % Get association stat
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    bhv(curr_animal_idx).move_rate = move_rate;
    bhv(curr_animal_idx).move_rate_stim = move_rate_stim;
    bhv(curr_animal_idx).move_rate_nonstim = move_rate_nonstim;
    bhv(curr_animal_idx).p_cue_given_move = p_cue_given_move;
    bhv(curr_animal_idx).rxn_stat_p = rxn_stat_p;

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});
    AP_print_progress_fraction(curr_recording,length(recordings));

end

ld = cellfun(@(x) ((1:size(x,1)) - find(x<0.05,1))',{bhv.rxn_stat_p}','uni',false);
use_animals = ~cellfun(@isempty,ld);

move_rate_cat = cell2mat({bhv(use_animals).move_rate}');
move_rate_stim_cat = cell2mat({bhv(use_animals).move_rate_stim}');
move_rate_nonstim_cat = cell2mat({bhv(use_animals).move_rate_nonstim}');
p_c2m_cat = cell2mat({bhv(use_animals).p_cue_given_move}');

ld_cat = cell2mat(ld(use_animals));
ld_split = ld_cat + linspace(0,(n_bins-1)/n_bins,n_bins);

[move_rate_avg,move_rate_grp] = ap.groupfun(@mean,move_rate_cat(:),ld_split(:));
move_rate_sem = ap.groupfun(@AP_sem,move_rate_cat(:),ld_split(:));

[move_rate_stim_avg,move_rate_stim_grp] = ap.groupfun(@mean,move_rate_stim_cat(:),ld_split(:));
move_rate_stim_sem = ap.groupfun(@AP_sem,move_rate_stim_cat(:),ld_split(:));

[move_rate_nonstim_avg,move_rate_nonstim_grp] = ap.groupfun(@mean,move_rate_nonstim_cat(:),ld_split(:));
move_rate_nonstim_sem = ap.groupfun(@AP_sem,move_rate_nonstim_cat(:),ld_split(:));

[p_c2m_avg,p_c2m_avg_grp] = ap.groupfun(@mean,p_c2m_cat(:),ld_split(:));
p_c2m_sem = ap.groupfun(@AP_sem,p_c2m_cat(:),ld_split(:));

figure('Name','Fig S1 move stim rate prob'); tiledlayout(2,1);
ax1 = nexttile; hold on;
h1 = errorbar( ...
    padarray(reshape(move_rate_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_sem,n_bins,[]),[1,0],NaN,'post'), ...
    'color',[0,0,0],'linewidth',2,'capsize',0);
h2 = errorbar( ...
    padarray(reshape(move_rate_stim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_sem,n_bins,[]),[1,0],NaN,'post'), ...
    'color',[0,0.5,0],'linewidth',2,'capsize',0);
h3 = errorbar( ...
    padarray(reshape(move_rate_nonstim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_sem,n_bins,[]),[1,0],NaN,'post'), ...
    'color',[0.5,0,0],'linewidth',2,'capsize',0);
xlabel('Learned day')
ylabel(' Wheel turns/s');
xline(0,'r');
legend([h1(1),h2(1),h3(1)],{'All','Stim','Non-stim'});

ax2 = nexttile;
errorbar( ...
    padarray(reshape(p_c2m_avg_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_sem,n_bins,[]),[1,0],NaN,'post'),'k','linewidth',2);
xlabel('Learned day')
ylabel('P(stim|move)');
xline(0,'r');

linkaxes([ax1,ax2],'x');
ap.prettyfig;
xlim([-3.1,2.9])

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end

% Specify re-load for next data (this fig uses raw instead of pre-packaged)
load_dataset_overwrite = true;


%% [Fig S1F] Fraction of quiescent passive trials

%%% Load non-activity data
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Get which trials in passive were quiescent
animals = unique(bhv.animal,'stable');

passive_quiescent_trials = struct( ...
    'stim_x',cell(length(animals),1), ...
    'quiescent',cell(length(animals),1));

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

    for use_rec = 1:length(ephys_days)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        load_parts.behavior = true;
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

        passive_quiescent_trials(animal_idx).stim_x{use_rec,1} = trial_stim_values;
        passive_quiescent_trials(animal_idx).quiescent{use_rec,1} = quiescent_trials;

    end
    ap.print_progress_fraction(animal_idx,length(animals));
end

stim_unique = unique(passive_quiescent_trials(1).stim_x{1});
stim_color = ap.colormap('BKR',length(stim_unique));

% Plot quiescent trials pre/post learning
passive_quiescent_stim = cell2mat(cellfun(@(stim,quiescent) ...
    ap.groupfun(@mean,+quiescent,stim)', ...
    cat(1,passive_quiescent_trials.stim_x), ...
    cat(1,passive_quiescent_trials.quiescent),'uni',false));

[passive_quiescent_stim_ld,passive_quiescent_stim_ld_grp] = ...
    ap.nestgroupfun({@mean,@mean},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);
passive_quiescent_stim_ld_sem = ...
    ap.nestgroupfun({@mean,@ap.sem},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);

figure('Name','Fig S1 frac q trials');
x_labels = ["Pre-learn","Post-learn"];
errorbar(reordercats(categorical(x_labels),x_labels), ...
    passive_quiescent_stim_ld,passive_quiescent_stim_ld_sem, ...
    'linewidth',2,'capsize',0)
axis padded;
ylabel('Frac. quiescent trials');
set(gca,'ColorOrder',stim_color)
ap.prettyfig;

% Plot quiescent trials binned by trial number
n_trialsplit = length(stim_unique)*5; % n presentations of each stim
passive_trialsplit_idx = cell2mat(cellfun(@(stim) ...
    (floor((0:length(stim)-1)/n_trialsplit)+1)', ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

passive_ld_idx = cell2mat(cellfun(@(rec,stim) repelem(rec,length(stim))', ...
    num2cell(bhv.days_from_learning), ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

[passive_quiescent_stim_trialsplit,passive_quiescent_stim_trialsplit_grp] = ...
    ap.groupfun(@mean,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_sem = ...
    ap.groupfun(@ap.sem,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit,[],[],nan);

passive_quiescent_stim_trialsplit_sem_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit_sem,[],[],nan);

passive_quiescent_stim_daysplit_x = reshape((unique(passive_quiescent_stim_trialsplit_grp(:,1)) + ...
    linspace(0,1,max(passive_quiescent_stim_trialsplit_grp(:,2)+1)))',[],1);

figure('Name','Fig S1 daysplit q trials');
plot_x_idx = isbetween(passive_quiescent_stim_daysplit_x,-3,3);

passive_quiescent_stim_trialsplit_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_grid,[0,1],nan,'post'),[2,1,3]),[],3);
passive_quiescent_stim_trialsplit_sem_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_sem_grid,[0,1],nan,'post'),[2,1,3]),[],3);

ap.errorfill(passive_quiescent_stim_daysplit_x(plot_x_idx), ...
    passive_quiescent_stim_trialsplit_grid_reshape(plot_x_idx,:), ...
    passive_quiescent_stim_trialsplit_sem_grid_reshape(plot_x_idx,:),stim_color);
xline(0,'k');
xlabel('Days from learning');
ylabel('Frac. quiescent trials');
axis padded;
ap.prettyfig;

% ~~~ STATS ~~~
print_stat('\n--FIG S1--\n');

sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

[passive_quiescent_stim_ld_animal,passive_quiescent_stim_ld_animal_grp] = ...
    ap.groupfun(@mean,passive_quiescent_stim, ...
    [grp2idx(bhv.animal),bhv.days_from_learning >= 0]);

for curr_stim = 1:length(stim_unique)
    stat_p = ranksum( ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 0,curr_stim), ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 1,curr_stim));

    print_stat('Ranksum pre/post stim %d: p = %.2g%s\n', ...
        stim_unique(curr_stim),stat_p,sig_flag(stat_p));
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S2A] Histology slices

animals = [ ...
    "AM011","AM012","AM014","AM015","AM016","AM017", ...
    "AM018","AM019","AM021","AM022","AM026","AM029", ...
    "AP023","AP025"];

% (manually selected slice with probe tract for each animal)
animal_slices = [ ...
    6, 5, 8, 2, 5, 8, ...
    6, 5, 3, 3, 4, 3, ...
    4, 3];

figure('Name','Fig S2 histology'); h = tiledlayout('flow','TileSpacing','tight');

for curr_animal_idx = 1:length(animals)

    animal = animals(curr_animal_idx);
    use_slice = animal_slices(curr_animal_idx);

    % Just load all images
    histology_path = plab.locations.filename('server',animal,[],[],'histology');
    histology_dir = dir(fullfile(histology_path,'*.tif'));

    histology_filenames = cellfun(@(path,name) fullfile(path,name), ...
        {histology_dir.folder},{histology_dir.name},'uni',false);
    [~,sort_idx] = natsortfiles(histology_filenames);

    histology_im = tiffreadVolume(histology_filenames{sort_idx(use_slice)});
   
    % Plot slice as RGB
    chan_cols = [0,1,0;1,0,0];

    chan_clim = double(prctile(histology_im,[10,90],[1,2])).*[1;2];

    histology_im_rgb = min(1,sum(cell2mat(arrayfun(@(chan) ...
        mat2gray(histology_im(:,:,chan),chan_clim(:,:,chan)).* ...
        permute(chan_cols(chan,:),[1,3,2]), ...
        permute(1:size(histology_im,3),[1,3,4,2]),'uni',false)),4));

    nexttile; 
    downsample_factor = 5;
    image(imresize(histology_im_rgb,1/downsample_factor)); 
    axis image off;

end

ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S2B] CCF-aligned probe histology

[av,tv,st] = ap_histology.load_ccf;

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Set atlas bins to plot through striatum
striatum_ccf_id = find(contains(lower(st.safe_name),'caudoputamen'));
striatum_ccf_ap = prctile(find(max(av == striatum_ccf_id,[],[2,3])),[0,100]);

n_atlas_bins = 10;
atlas_bins = round(linspace(striatum_ccf_ap(1),striatum_ccf_ap(2),n_atlas_bins+1));

histology_atlas_bin_max = cell(length(animals),n_atlas_bins);

for curr_animal_idx = 1:length(animals)

    histology_path = plab.locations.filename('server',animals{curr_animal_idx},[],[],'histology');
    load(fullfile(histology_path,'AP_histology_processing.mat'));

    % Load images
    image_path = histology_path;
    image_dir = dir(fullfile(image_path,'*.tif'));
    image_filenames = cellfun(@(path,name) fullfile(path,name), ...
        {image_dir.folder},{image_dir.name},'uni',false);
    [~,sort_idx] = ap_histology.natsortfiles(image_filenames);

    images = cell(size(image_dir));
    for curr_im = 1:length(sort_idx)
        images{curr_im} = tiffreadVolume( ...
            image_filenames{sort_idx(curr_im)});
    end

    % Grab atlas images
    n_slices = length(images);
    slice_atlas = struct('tv',cell(n_slices,1), 'av',cell(n_slices,1));
    slice_atlas_ccf = struct('ap',cell(n_slices,1),'ml',cell(n_slices,1),'dv',cell(n_slices,1));
    for curr_slice = 1:length(images)
        [slice_atlas(curr_slice),slice_atlas_ccf(curr_slice)] = ...
            ap_histology.grab_atlas_slice(av,tv, ...
            AP_histology_processing.histology_ccf.slice_vector, ...
            AP_histology_processing.histology_ccf.slice_points(curr_slice,:), 1);
    end

    % Build volume of histology images
    histology_volume = zeros(size(tv),'single');
    probe_channel = 2;
    for curr_im_idx = 1:length(images)

        % Rigid transform
        im_rigid_transformed = ap_histology.rigid_transform( ...
            images{curr_im_idx}(:,:,probe_channel),curr_im_idx,AP_histology_processing);

        % Affine/nonlin transform
        if isfield(AP_histology_processing.histology_ccf,'control_points') && ...
                (size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) == ...
                size(AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},1)) && ...
                size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) >= 3
            % Manual alignment (if >3 matched points)
            histology2atlas_tform = fitgeotform2d( ...
                AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx}, ...
                AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},'pwl');
        elseif isfield(AP_histology_processing.histology_ccf,'atlas2histology_tform')
            % Automatic alignment
            histology2atlas_tform = invert(AP_histology_processing.histology_ccf.atlas2histology_tform{curr_im_idx});
        end

        atlas_slice_aligned = imwarp(im_rigid_transformed, ...
            histology2atlas_tform,'nearest','OutputView', ...
            imref2d(size(slice_atlas(curr_im_idx).av)));

        % % Check match
        % figure; imshowpair(slice_atlas(curr_im_idx).av,atlas_slice_aligned);

        % Add points to volume in CCF space
        curr_ccf_idx = sub2ind(size(tv), ...
            round(slice_atlas_ccf(curr_im_idx).ap(:)), ...
            round(slice_atlas_ccf(curr_im_idx).dv(:)), ...
            round(slice_atlas_ccf(curr_im_idx).ml(:)));

        histology_volume(curr_ccf_idx) = histology_volume(curr_ccf_idx) + ...
            single(atlas_slice_aligned(:));

    end

    % Get max of histology volume in atlas bins
    for curr_atlas_bin = 1:n_atlas_bins
        plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
        histology_atlas_bin_max{curr_animal_idx,curr_atlas_bin} = ...
            permute(max(histology_volume(atlas_bins(curr_atlas_bin): ...
            atlas_bins(curr_atlas_bin+1),:,:),[],1),[2,3,1]);
    end
    
    ap.print_progress_fraction(curr_animal_idx,length(animals));
end

% Plot probe channel overlay by colored animals
animal_colors = ap.colormap('tube',14);
overlay_dilation = 1;
histology_clim = repelem({[200,500]},length(animals),1);

figure('Name','Fig S2 colored histology'); tiledlayout('TileSpacing','none');
for curr_atlas_bin = 1:n_atlas_bins

    % Plot CCF borders from the middle of the bin
    plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
    curr_ccf_borders = imdilate(boundarymask(permute(av(plot_atlas_ap,:,:),[2,3,1])),ones(overlay_dilation));

    % Max, color, and flip contrast
    curr_histology_volume_max_gray = cellfun(@(x,c) ...
        mat2gray(x,c),histology_atlas_bin_max(:,curr_atlas_bin),histology_clim,'uni',false);
        
    curr_histology_volume_max_colored = ...
        cat(4,curr_histology_volume_max_gray{:}).*permute(animal_colors,[3,4,2,1]);
    
    curr_histology_volume_max_colored_white = ...
        curr_histology_volume_max_colored+(1-cat(4,curr_histology_volume_max_gray{:}));

    curr_histology_combined = min(curr_histology_volume_max_colored_white,[],4);

    % Plot CCF over combined colored image
    curr_overlay = imoverlay(curr_histology_combined,curr_ccf_borders,'k');
    nexttile;imagesc(curr_overlay);axis image off;
    drawnow;
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S3] Striatal domain clustering and classification

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

% Choose animal and day to plot
use_animal = 'AM026';
use_ld = 0;
use_rec = strcmp(bhv.animal,use_animal) & bhv.days_from_learning == use_ld;
use_cortex_kernel = ctx_str_maps.cortex_striatum_map{use_rec};

domain_color = {'180','60','310'};

% Plot all domains
% (grayscale and colored by domain)
figure('Name','Fig S3 map examples');
h = tiledlayout(size(use_cortex_kernel,3),2,'TileSpacing','none');
for curr_depth=1:size(use_cortex_kernel, 3)
    nexttile;
    imagesc(use_cortex_kernel(:,:,curr_depth));
    axis image off;
    ap.wf_draw('cortex',[0.5,0.5,0.5]);
    colormap(gca,ap.colormap('WK',[],2));

    nexttile;
    imagesc(use_cortex_kernel(:,:,curr_depth));
    axis image off
    ap.wf_draw('cortex',[0.5,0.5,0.5]);
    colormap(gca,ap.colormap(['W' domain_color{domain_idx_rec{use_rec}(curr_depth)}],[],2));
end
clim(h.Children,[0,0.01]);
ap.prettyfig;

% Plot domain means and ROIs
striatum_wf_roi_outlines = ...
    arrayfun(@(curr_domain) ...
    bwboundaries(striatum_wf_roi(:,:,curr_domain)),1:n_domains);

figure('Name','Fig S3 domain maps');
h = tiledlayout(n_domains,1,'tilespacing','none');
for curr_domain = 1:n_domains
    nexttile;
    imagesc(kmeans_cluster_mean(:,:,curr_domain));
    axis image off
    colormap(gca,ap.colormap(['W',domain_color{curr_domain}],[],2));
    ap.wf_draw('cortex',[0.5,0.5,0.5]);

    patch(striatum_wf_roi_outlines{curr_domain}(:,2), ...
        striatum_wf_roi_outlines{curr_domain}(:,1),'w', ...
        'FaceColor','none','EdgeColor','y','linewidth',2);
end
clim(h.Children,[0,0.01]);
ap.prettyfig;

% Load data from example recording
animal = bhv.animal{use_rec};
rec_day = bhv.rec_day{use_rec};
recordings = plab.find_recordings(animal,rec_day,'stim_wheel*');
rec_time = recordings.recording{end};
load_parts.ephys = true;
ap.load_recording;

% Plot units and overlay clustering
figure('Name','Fig S3 example units'); 
ax = axes; hold on;

domain_color_rgb = [1,0,0;0,1,0;0,0,1];
domain_im = permute(domain_color_rgb(domain_idx_rec{use_rec},:),[1,3,2]);
imagesc(ax,[],ctx_str_maps.depth_group_edges{use_rec},domain_im);
ax.YDir = 'reverse';

ap.plot_unit_depthrate(spike_times_timelite,spike_templates,template_depths,[],ax);
yline(ctx_str_maps.depth_group_edges{use_rec},'linewidth',2,'color',[0.5,0.5,0.5]);

ap.prettyfig;

% Plot clustered MUA
depth_group = discretize(spike_depths,ctx_str_maps.depth_group_edges{use_rec});

plot_t = [100,140];
bin_t = 0.1;
domain_mua = zeros(n_domains,diff(plot_t)/bin_t);
for curr_domain = 1:n_domains
    curr_spikes = spike_times_timelite(isbetween(spike_times_timelite,plot_t(1),plot_t(2)) & ...
        ismember(depth_group,find(domain_idx_rec{use_rec} == curr_domain)));
    domain_mua(curr_domain,:) = histcounts(curr_spikes,plot_t(1):bin_t:plot_t(2))./bin_t;
end

figure('Name','Fig S3 example mua');
plot(conv(plot_t(1):bin_t:plot_t(2),[0.5,0.5],'valid'),domain_mua','linewidth',2);
set(gca,'ColorOrder',domain_color_rgb);
axis off;
ap.scalebar(10,500);
ap.prettyfig

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S4] Task activity heatmaps (Fig 1A) split by animal

%%% Load non-activity data
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

% plot_day_bins = [-Inf,-2,0,Inf];
plot_day_bins = [-Inf,-3:2,Inf];

cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Plot single area by animal
animals = unique(bhv.animal,'stable');

plot_domains = [1,2];

trials_scale = [15,50]; % (small, big)

heatmap_smooth = [1,1]; % (don't smooth for animal split)

for curr_domain = plot_domains

    f_ctx = figure('Name',sprintf('Fig S4 cortex %d heatmaps',curr_domain));
    h_ctx = tiledlayout(f_ctx,length(animals),max(cortex_plot_day_grp), ...
        'TileIndexing','ColumnMajor','TileSpacing','none');

    f_str = figure('Name',sprintf('Fig S4 striatum %d',curr_domain));
    h_str = tiledlayout(f_str,length(animals),max(cortex_plot_day_grp), ...
        'TileIndexing','ColumnMajor','TileSpacing','none');

    for curr_day = 1:max(cortex_plot_day_grp)
        for curr_animal = 1:length(animals)

            % Cortex
            curr_trials = find(cortex_plot_day_grp == curr_day & ...
                wf_grp.animal == curr_animal);

            if any(curr_trials)
                nexttile(h_ctx,tilenum(h_ctx,curr_animal,curr_day));

                [sorted_rxn,sort_idx] = sort(wf_grp.rxn(curr_trials));
                imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_domain,1),heatmap_smooth));
                colormap(gca,ap.colormap('WG'));
                clim([0,1e-2]);
                axis off;

                hold on
                xline(0,'color','r');
                plot(sorted_rxn,1:length(curr_trials),'b');

                h = ap.scalebar([],trials_scale(find(length(curr_trials) > trials_scale,1,'last')));
                if length(curr_trials) < trials_scale(2)
                    h(2).Color = [0,1,1];
                end
            end

            % Striatum
            curr_trials = find(striatum_plot_day_grp == curr_day & ...
                striatum_mua_grp.domain_idx == curr_domain & ...
                striatum_mua_grp.animal == curr_animal);

            if any(curr_trials)
                nexttile(h_str,tilenum(h_str,curr_animal,curr_day));

                [sorted_rxn,sort_idx] = sort(striatum_mua_grp.rxn(curr_trials));
                imagesc(psth_t,[],movmean(striatum_mua(curr_trials(sort_idx),:,1),heatmap_smooth));
                colormap(gca,ap.colormap('WK'));
                clim([0,3]);
                axis off;

                hold on
                xline(0,'color','r');
                plot(sorted_rxn,1:length(curr_trials),'b');

                h = ap.scalebar([],trials_scale(find(length(curr_trials) > trials_scale,1,'last')));
                if length(curr_trials) < trials_scale(2)
                    h(2).Color = [0,1,1];
                end
            end

            drawnow;
        end
    end
    ap.prettyfig('eps',f_ctx);
    ap.prettyfig('eps',f_str);
end

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S5] Non-stim move activity

%%% Load non-activity data
load_dataset = 'task';
Marica_2026.figures.load_data;
%%%

% Load nonstim move activity
U_master = plab.wf.load_master_U;
load(fullfile(data_path,'nonstim_move'));

use_nostim_move_recordings = ~cellfun(@isempty,nonstim_move.V_move_nostim_align);
nostim_move_wheel_t = nonstim_move.wheel_align_time{find(use_nostim_move_recordings,1)};

%%% NON-STIM MOVE DATA PREPROCESSING
% Get widefield ROIs for no stim moves
wf_nostim_move_striatum_roi = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master,x,[],[],striatum_wf_roi), ...
    [3,2,1]),nonstim_move.V_move_nostim_align(use_nostim_move_recordings), ...
    'uni',false));
baseline_t = wf_t < 0;
wf_nostim_move_striatum_roi = wf_nostim_move_striatum_roi - nanmean(wf_nostim_move_striatum_roi(:,baseline_t,:,:),2);

% Get nonstim move ephys
% (sum into domain multiunit)
striatum_nostim_move_mua_sum = cellfun(@(mua,domain_idx) ...
    permute(ap.groupfun(@sum,mua,domain_idx,[]),[3,2,1]), ...
    nonstim_move.binned_msn_spikes_move_nostim_align,domain_idx_rec,'uni',false);

% (smooth and normalize) 
baseline_t = psth_t < -0.3;
striatum_nostim_baseline = cellfun(@(mua) ...
        mean(mua(:,baseline_t,:,1),[2]), ...
        striatum_nostim_move_mua_sum,'uni',false, ...
        'ErrorHandler',@(varargin) NaN);
striatum_nostim_move_mua = cellfun(spikes_norm_smooth_reshape_fcn, ...
        striatum_nostim_move_mua_sum,striatum_nostim_baseline,'uni',false); %#ok<FUNFUN>

% Plot move-aligned (stim, non-stim) binned by day
% (non-stim not saved by trial, so day bins by weighted average)
plot_day_bins = [-Inf,0,Inf];

plot_day_colors = [0,0,0;0.7,0,0];
plot_day_color_letters = {'k','r'}; % for ap.colormap in rxn plot

day_grp = discretize(max(bhv.days_from_learning,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

figure('Name','Fig S5 nostim striatum'); h_striatum = tiledlayout(n_domains,3); title(h_striatum,'Striatum');
figure('Name','Fig S5 nostim cortex'); h_cortex = tiledlayout(n_domains,3); title(h_cortex,'Cortex');
figure('Name','Fig S5 nostim rxn'); h_rxn = tiledlayout(2,1);

figure('Name','Fig S5 nostim wheel'); h_wheel = tiledlayout(1,3);
for curr_plot = 1:3
    nexttile(h_wheel,curr_plot); hold on; set(gca,'ColorOrder',plot_day_colors);
end

for curr_domain = 1:n_domains

    h_striatum_sub = gobjects(3,1);
    h_cortex_sub = gobjects(3,1);
    for sub = 1:3
        h_striatum_sub(sub) = nexttile(h_striatum); hold on;
        set(gca,'ColorOrder',plot_day_colors);

        h_cortex_sub(sub) = nexttile(h_cortex); hold on;
        set(gca,'ColorOrder',plot_day_colors);
    end

    for curr_day_grp = 1:length(plot_day_bins)-1

        % Get day-binned (weighted average) activity
        curr_nomove_idx = use_nostim_move_recordings & day_grp==curr_day_grp;
        n_moves = cellfun(@(x) size(x,1),nonstim_move.move_nostim_wheel(curr_nomove_idx));

        animals = unique(bhv.animal);
        [~,curr_nostim_move_animal_idx] = ismember(nonstim_move.animal(curr_nomove_idx),animals);

        % (striatum)
        curr_nostim_striatum = cellfun(@(x,domain) x(domain==curr_domain,:), ...
            striatum_nostim_move_mua(curr_nomove_idx), ...
            striatum_domain_idx(curr_nomove_idx),'uni',false);
        
        curr_striatum_nostim_weighted = cellfun(@(act,n) act.*n,curr_nostim_striatum,num2cell(n_moves),'uni',false);
        
        curr_striatum_nostim_wavg = ap.groupfun(@sum,cell2mat(curr_striatum_nostim_weighted), ...
            curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)))./ ...
            ap.groupfun(@sum,n_moves(~cellfun(@isempty,curr_striatum_nostim_weighted)), ...
             curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)));

        % (cortex)
        curr_nostim_cortex = wf_nostim_move_striatum_roi( ...
            day_grp(use_nostim_move_recordings)==curr_day_grp,:,curr_domain);

        curr_nostim_cortex_wavg = ...
            ap.groupfun(@sum,curr_nostim_cortex.*n_moves,curr_nostim_move_animal_idx)./ ...
            ap.groupfun(@sum,n_moves,curr_nostim_move_animal_idx);
        
        % Get stim move activity 
        % (only use animals from above to pair data)
        
        % (striatum)
        curr_nomove_animals_striatum = unique(curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)));
        curr_striatum_trials_idx = ismember(striatum_mua_grp.animal,curr_nomove_animals_striatum) & ...
            striatum_mua_grp.domain_idx == curr_domain & ...
            striatum_plot_day_grp == curr_day_grp;

        curr_striatum_stimmove = ap.groupfun(@mean, ...
            striatum_mua(curr_striatum_trials_idx,:,2), ...
            striatum_mua_grp.animal(curr_striatum_trials_idx));

        % (cortex) 
        curr_nomove_animals_cortex = unique(curr_nostim_move_animal_idx);
        curr_cortex_trials_idx = ismember(wf_grp.animal,curr_nomove_animals_cortex) & ...
            cortex_plot_day_grp == curr_day_grp;
       
        curr_cortex_stimmove = ap.groupfun(@mean, ...
            wf_striatum_roi(curr_cortex_trials_idx,:,curr_domain,2), ...
            wf_grp.animal(curr_cortex_trials_idx));

        % Plot
        axes(h_striatum_sub(1));ap.errorfill(psth_t,nanmean(curr_striatum_stimmove,1),ap.sem(curr_striatum_stimmove,1));
        axes(h_striatum_sub(2));ap.errorfill(psth_t,nanmean(curr_striatum_nostim_wavg,1),ap.sem(curr_striatum_nostim_wavg,1));
        axes(h_striatum_sub(3));ap.errorfill(psth_t,nanmean(curr_striatum_stimmove-curr_striatum_nostim_wavg,1), ...
            ap.sem(curr_striatum_stimmove-curr_striatum_nostim_wavg,1));

        axes(h_cortex_sub(1));ap.errorfill(wf_t,nanmean(curr_cortex_stimmove,1),ap.sem(curr_cortex_stimmove,1));
        axes(h_cortex_sub(2));ap.errorfill(wf_t,nanmean(curr_nostim_cortex_wavg,1),ap.sem(curr_nostim_cortex_wavg,1));
        axes(h_cortex_sub(3));ap.errorfill(wf_t,nanmean(curr_cortex_stimmove-curr_nostim_cortex_wavg,1), ...
            ap.sem(curr_cortex_stimmove-curr_nostim_cortex_wavg,1));

        if curr_domain == 1
            
            % Plot histogram of stim relative to move onset (on first domain)
            rxn_bins = [-0.5:0.01:1];
            rxn_bin_x = [rxn_bins(2),rxn_bins(2:end-2)+diff(rxn_bins(2:end-1))/2,rxn_bins(end-1)];
            rxn_histogram = cell2mat(arrayfun(@(x) histcounts(-striatum_mua_grp.rxn(curr_striatum_trials_idx & ...
                striatum_mua_grp.animal == x),rxn_bins,'Normalization','probability'), ...
                (1:length(unique(striatum_mua_grp.animal)))','uni',false));
            nexttile(h_rxn,curr_day_grp);
            imagesc(rxn_bin_x,[],nanmean(rxn_histogram,1));
            colormap(nexttile(h_rxn,curr_day_grp),ap.colormap(...
                sprintf('W%s',upper(plot_day_color_letters{curr_day_grp}))));
            set(nexttile(h_rxn,curr_day_grp),'XTick','','YTick','');
            box on;colorbar;

            % Plot wheel
            wheel_stim_animal = arrayfun(@(x) nanmean(cell2mat(nonstim_move.move_stim_wheel( ...
                (grp2idx(bhv.animal) == x) & (day_grp == curr_day_grp))),1), ...
                unique(grp2idx(bhv.animal)),'uni',false);  
            nexttile(h_wheel,1); title('Stim');
            ap.errorfill(nonstim_move.wheel_align_time{end}, ...
                nanmean(vertcat(wheel_stim_animal{:}),1),ap.sem(vertcat(wheel_stim_animal{:}),1));    

            wheel_nostim_animal = arrayfun(@(x) nanmean(cell2mat(nonstim_move.move_nostim_wheel( ...
                (grp2idx(bhv.animal) == x) & (day_grp == curr_day_grp))),1), ...
                unique(grp2idx(bhv.animal)),'uni',false);        
            nexttile(h_wheel,2); title('No stim');
            ap.errorfill(nonstim_move.wheel_align_time{end}, ...
                nanmean(vertcat(wheel_nostim_animal{:}),1),ap.sem(vertcat(wheel_nostim_animal{:}),1));  

            nexttile(h_wheel,3); title('Difference');
            wheel_diff_animal = cellfun(@(stim_wheel,nostim_wheel) stim_wheel-nostim_wheel, ...
                wheel_stim_animal,wheel_nostim_animal,'uni',false);
            ap.errorfill(nonstim_move.wheel_align_time{end}, ...
                nanmean(vertcat(wheel_diff_animal{:}),1),ap.sem(vertcat(wheel_diff_animal{:}),1));    

            % Draw activity scalebars
            axes(h_cortex_sub(1)); ap.scalebar(0.5,5e-3);
            axes(h_striatum_sub(1)); ap.scalebar(0.5,1);
            axes(nexttile(h_wheel,1)); deg_scale = 360; ap.scalebar(0.5,deg_scale*(1024/360)/3);

        end

    end
end
% (link x/y axes for activity)
linkaxes(h_striatum.Children,'xy');
linkaxes(h_cortex.Children,'xy');
linkaxes(h_wheel.Children,'xy');

% (draw lines at t = 0)
arrayfun(@(x) xline(x,0),[h_striatum.Children])
arrayfun(@(x) xline(x,0),[h_cortex.Children])
arrayfun(@(x) xline(x,0),[h_wheel.Children])
arrayfun(@(x) xline(x,0),[h_rxn.Children(end:-2:1)])

% (match color axes for stim histogram)
clim(nexttile(h_rxn,1),clim(nexttile(h_rxn,2)));

ap.prettyfig([],h_striatum.Parent);
ap.prettyfig([],h_cortex.Parent);
ap.prettyfig([],h_wheel.Parent);
ap.prettyfig([],h_rxn.Parent);

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% [Fig S9] Striatum cell type properties 

% Load ephys properties
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2026','data');
load(fullfile(data_path,'ephys_properties'));

% Concatenate data
striatum_celltypes = ["msn","fsi","tan"];

striatum_celltype_cat.msn = logical(vertcat(ephys_properties.str_msn_idx{:}));
striatum_celltype_cat.fsi = logical(vertcat(ephys_properties.str_fsi_idx{:}));
striatum_celltype_cat.tan = logical(vertcat(ephys_properties.str_tan_idx{:}));

waveform_norm_cat = vertcat(ephys_properties.waveform{:})./max(abs(vertcat(ephys_properties.waveform{:})),[],2);

acg_cat = vertcat(ephys_properties.acg{:});
acg_t = -500:500; % (just hardcoding - it's set somewhere in bombcell)

waveform_duration_cat = vertcat(ephys_properties.waveformDuration_peakTrough_us{:});
postspike_suppression_cat = vertcat(ephys_properties.postSpikeSuppression_ms{:});
firing_rate_cat = vertcat(ephys_properties.mean_firingRate{:});

% (remove units that are likely light artifacts)
framerate = 70;
light_artifact_units = max(acg_cat(:,acg_t<-200),[],2) > framerate;
for curr_celltype = striatum_celltypes
    striatum_celltype_cat.(curr_celltype)(light_artifact_units) = false;
end

figure('Name','Fig S9 celltype features');
h = tiledlayout(length(striatum_celltypes),2);
for curr_celltype = striatum_celltypes
    nexttile;
    ap.errorfill([],nanmean(waveform_norm_cat(striatum_celltype_cat.(curr_celltype),:),1), ...
        nanstd(waveform_norm_cat(striatum_celltype_cat.(curr_celltype),:),1));

    nexttile;
    ap.errorfill(acg_t,nanmean(acg_cat(striatum_celltype_cat.(curr_celltype),:),1), ...
        nanstd(acg_cat(striatum_celltype_cat.(curr_celltype),:),1));
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy');
xlim(h.Children(1),[-1,1].*300);
ylim(h.Children(1),[0,50]);
ap.prettyfig;

% Histograms of properties
figure('Name','Fig S9 celltype histograms'); tiledlayout(length(striatum_celltypes),3);
for curr_celltype = striatum_celltypes
    nexttile; 
    histogram(waveform_duration_cat(striatum_celltype_cat.(curr_celltype)), ...
        [-Inf,linspace(0,prctile(waveform_duration_cat,99),20),Inf], ...
        'normalization','probability', ...
        'FaceColor','k','FaceAlpha',1,'EdgeColor','none');
    xlabel('Waveform duration'); xline(400,'r','linewidth',2);
    title(curr_celltype);

    nexttile; 
    histogram(postspike_suppression_cat(striatum_celltype_cat.(curr_celltype)), ...
        [-Inf,linspace(0,prctile(postspike_suppression_cat,99),20),Inf], ...
        'normalization','probability', ...
        'FaceColor','k','FaceAlpha',1,'EdgeColor','none');
    xlabel('Postspike suppression'); xline(40,'r','linewidth',2);
    title(curr_celltype);

    nexttile; 
    histogram(firing_rate_cat(striatum_celltype_cat.(curr_celltype)), ...
        [-Inf,linspace(0,prctile(firing_rate_cat,99),20),Inf], ...
        'normalization','probability', ...
        'FaceColor','k','FaceAlpha',1,'EdgeColor','none');
    xlabel('Firing rate');
    title(curr_celltype);
end
ap.prettyfig;

% Plot number/fraction of units
celltype_rec_n = cell2mat(arrayfun(@(x) cellfun(@(celltype,str) sum(celltype(str)), ...
    ephys_properties.(sprintf('str_%s_idx',x)), ...
    ephys_properties.striatal_units),striatum_celltypes,'uni',false));

celltype_rec_frac = cell2mat(arrayfun(@(x) cellfun(@(celltype,str) mean(celltype(str)), ...
    ephys_properties.(sprintf('str_%s_idx',x)), ...
    ephys_properties.striatal_units),striatum_celltypes,'uni',false));

figure('Name','Fig S9 celltype n'); tiledlayout(1,2);
nexttile; hold on;
bar(striatum_celltypes,nanmean(celltype_rec_n,1));
errorbar(categorical(striatum_celltypes),nanmean(celltype_rec_n,1), ...
    nanstd(celltype_rec_n,1),'k', ...
    'LineStyle','none','LineWidth',2,'CapSize',0);
ylabel('Number of units')

nexttile; hold on;
bar(striatum_celltypes,nanmean(celltype_rec_frac,1));
errorbar(categorical(striatum_celltypes),nanmean(celltype_rec_frac,1), ...
    nanstd(celltype_rec_frac,1),'k', ...
    'LineStyle','none','LineWidth',2,'CapSize',0);
ylabel('Fraction of units');

ap.prettyfig;

% ~~~ SAVE FIGS ~~~
if exist('fig_save_flag','var') && fig_save_flag
    save_figs();
    close(findall(0,'Type','figure'));
end


%% Close and clear

% Close stat file
if exist('stat_fid','var')
    fclose(stat_fid);
end

clearvars

toc







