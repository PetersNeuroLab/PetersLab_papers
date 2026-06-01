
%% Fig SB A-D: activity/encoding regressor examples
main_preload_vars = who;

animal = 'DS007';
rec_day = '2024-07-11';
rec_time = '0815';

load_parts.mousecam = true;
load_parts.widefield = true;
load_parts.widefield_master = true;
verbose = true;
ap.load_recording;

% Load and prep pupil SLEAP tracking
pupil_sleap_dir = dir(fullfile(fileparts(mousecam_fn),'**','*.h5'));
pupil_sleap_fn = fullfile(pupil_sleap_dir.folder,pupil_sleap_dir.name);

tracks = h5read(pupil_sleap_fn,'/tracks');       % frames x nodes x 2 (x/y-position)
% (scores aren't used at the moment - PG didn't need previously)
% pointScores = h5read(pupil_sleap_fn,'/point_scores')'; % frames x nodes
% instanceScores = h5read(pupil_sleap_fn, '/instance_scores')'; % transpose to 1 x frames

% Fit circle (solve for a,b,c in x^2 + y^2 + a*x + b*y + c = 0)
% (use only non-NaN vertices, remove points with too few vertices)
circle_fit_fun = @(x,y) [x(~any(isnan([x,y]),2)) y(~any(isnan([x,y]),2)) ...
    ones(sum(~any(isnan([x,y]),2)),1)]\ ...
    -(x(~any(isnan([x,y]),2)).^2+y(~any(isnan([x,y]),2)).^2);

min_pupil_points = 3;
pupil_valid_frames = sum(~any(isnan(tracks),3),2) >= min_pupil_points;

pupil_circle_fit = nan(3,length(mousecam_times));
pupil_circle_fit(:,pupil_valid_frames) = cell2mat(arrayfun(@(frame) ...
    circle_fit_fun(tracks(frame,:,1)',tracks(frame,:,2)'), ...
    find(pupil_valid_frames)','uni',false));

pupil = struct( ...
    'x',-pupil_circle_fit(1,:)./2, ...
    'y',-pupil_circle_fit(2,:)./2, ...
    'diameter',2*sqrt(sum(pupil_circle_fit(1:2,:).^2,1)/4-pupil_circle_fit(3,:)));

% Regression parameters
time_bins = [wf_t;wf_t(end)+1/wf_framerate];
n_components = 200;
frame_shifts = -10:30;
lamda_encode = 0;
lambda_decode = 15;
cv_fold = 5;
skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);

% Stim onset regressors
stim_regressors = histcounts(stimOn_times,time_bins);

% Move regressors
% (move onset, wheel velocity, pupil diameter, pupil velocity)
wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
move_onset_regressors = histcounts(wheel_starts,time_bins);

wheel_velocity_resample = interp1(timelite.timestamps,wheel_velocity,wf_t);
pupil_diameter_resample = interp1(mousecam_times,pupil.diameter,wf_t);
pupil_velocity_resample = interp1(sqrt(diff(pupil.x).^2+diff(pupil.y).^2),wf_t);

move_regressors = vertcat( ...
    move_onset_regressors,wheel_velocity_resample', ...
    pupil_diameter_resample',pupil_velocity_resample');

move_regressor_labels = {'Move onset','Wheel velocity', ...
    'Pupil diameter','Pupil velocity'};

% Regress movement or stimulus
regression_frames = skip_frames:size(wf_V,2)-skip_frames;

[kernels_move_encode,predicted_signals_move] = ...
    ap.regresskernel(move_regressors(:,regression_frames), ...
    wf_V(1:n_components,regression_frames), ...
    frame_shifts,lamda_encode,[],cv_fold);

[kernels_stim_encode,predicted_signals_stim] = ...
    ap.regresskernel(stim_regressors(:,regression_frames), ...
    wf_V(1:n_components,regression_frames), ...
    frame_shifts,lamda_encode,[],cv_fold);

% Plot raw and residual ROI activity
roi_mask_filename = fullfile(plab.locations.server_path,'Lab','Papers','Song_2025','data','General_information','roi.mat');
load(roi_mask_filename); % saved as `roi1`

plot_rois = {'l-V1','l-SSp'};
plot_roi_idx = ismember({roi1.name},plot_rois);
roi_masks = cell2mat(reshape(cellfun(@(x) x.mask,{roi1(plot_roi_idx).data},'uni',false),1,1,[]));
roi_labels = {roi1(plot_roi_idx).name};

wf_roi = ap.wf_roi(wf_U(:,:,1:n_components),wf_V(1:n_components,regression_frames),[],[],roi_masks);
wf_roi_moveresidual = ap.wf_roi(wf_U(:,:,1:n_components),wf_V(1:n_components,regression_frames)- ...
    fillmissing(predicted_signals_move,'constant',0),[],[],roi_masks);
wf_roi_stimresidual = ap.wf_roi(wf_U(:,:,1:n_components),wf_V(1:n_components,regression_frames)- ...
    fillmissing(predicted_signals_stim,'constant',0),[],[],roi_masks);

figure;
h = tiledlayout(size(move_regressors,1)+2,1,'tilespacing','none');
for curr_regressor = 1:size(move_regressors,1)
    nexttile;
    plot(wf_t(regression_frames),move_regressors(curr_regressor,regression_frames),'k','linewidth',2);
    ylabel(move_regressor_labels{curr_regressor});
end
for curr_roi = 1:size(roi_masks,3)
    nexttile;hold on;
    plot(wf_t(regression_frames),wf_roi(curr_roi,:),'color',[0,0.7,0],'linewidth',2);
    plot(wf_t(regression_frames),wf_roi_moveresidual(curr_roi,:),'k')
    plot(wf_t(regression_frames),wf_roi_stimresidual(curr_roi,:),'b')
    ylabel(roi_labels{curr_roi});
    legend({'Measured','Move-residual','Stim-residual'});
end
linkaxes(h.Children,'x');
xlim([870,895]);
ap.prettyfig;


%% Fig SB D-E: encoding and encoding-residual decoding

% Load master U
wf_U = plab.wf.load_master_U;
n_components = 200;

% Load behavior and kernels
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Song_2025','data');
load(fullfile(data_path,'behavior.mat'));
load(fullfile(data_path,'encoding_decoding_kernels.mat'));

% Concatenate workflows/learning days
workflow_cat = arrayfun(@(grp) behavior_each_mice{grp}.workflow_name,1:length(behavior_each_mice),'uni',false);
learned_cat = arrayfun(@(grp) behavior_each_mice{grp}.learned,1:length(behavior_each_mice),'uni',false);

% Set groups to use (VA,AV)
use_animal_grp = 1:2;
modality_idx_cat = cell2mat(cellfun(@(x) contains(x(~contains(x,'mixed')),'audio'), ...
    vertcat(workflow_cat{use_animal_grp}),'uni',false)); % 0 = Vis, 1 = Aud
learned_idx_cat = cell2mat(cellfun(@(x,workflow) x(~contains(workflow,'mixed')), ...
    vertcat(learned_cat{use_animal_grp}),vertcat(workflow_cat{use_animal_grp}),'uni',false));

% Get kernel types
kernel_types = string(fieldnames(encoding_decoding_kernels));

for use_kernel = kernel_types'
    % Grab kernel and average by modality/learning
    curr_kernels_animalsplit = [encoding_decoding_kernels(use_animal_grp).(use_kernel)];
    curr_kernels = horzcat(curr_kernels_animalsplit{:});

    use_kernels = ~cellfun(@isempty,curr_kernels);

    [kernel_avg,kernel_avg_grp] = ap.groupfun(@nanmean,cat(4, ...
        curr_kernels{use_kernels}),[],[],[], ...
        [modality_idx_cat(use_kernels),learned_idx_cat(use_kernels)]);

    % Get kernel pixels time-max
    kernels_px = plab.wf.svd2px(wf_U(:,:,1:n_components),kernel_avg);

     max_t = [-inf,inf]; % max over full kernel

    surround_samplerate = 35;
    frame_shifts = -10:30;
    frame_shifts_t = frame_shifts./surround_samplerate;
     % use_t = isbetween(frame_shifts_t,max_t(1),max_t(2)); 
     use_t = ones(size(frame_shifts_t));

    kernel_tmax = permute(max(kernels_px(:,:,use_t,:,:),[],3),[1,2,4,5,3]);

    % Plot time-max kernels
    plot_group_order = [0,0;0,1;1,0;1,1];
    [~,plot_grp_sort] = ismember(plot_group_order,kernel_avg_grp,'rows');
    plot_grp_order_name = {'Vis pre','Vis post','Aud pre','Aud post'};

    figure;
    h = tiledlayout(size(kernel_tmax,3),size(kernel_tmax,4),'TileSpacing','tight');
    modality_colors = {'WB','WR'};
    for curr_modal_learn = 1:size(kernel_tmax,4)
        for curr_subkernel = 1:size(kernel_tmax,3)
            % Choose tile
            curr_ax = nexttile(tilenum(h,curr_subkernel,curr_modal_learn));

            % Plot kernel tmax (in set group order)
            curr_plot_modal_learn = plot_grp_sort(curr_modal_learn);
            imagesc(kernel_tmax(:,:,curr_subkernel,curr_plot_modal_learn));
            axis image off;
            ap.wf_draw('ccf',[0.5,0.5,0.5]);

            curr_color = modality_colors{kernel_avg_grp(curr_plot_modal_learn,1)+1};
            colormap(curr_ax,ap.colormap(curr_color));

            if contains(use_kernel,'encode')
                clim([0,max(clim)]);
            elseif contains(use_kernel,'decode')
                clim([0,3e-4]);
            end
            colorbar

            % Title column
            if curr_subkernel == 1
                title(gca,plot_grp_order_name{curr_modal_learn})
            end      
        end
    end
    title(h,strrep(use_kernel,'_',' '));
    ap.prettyfig;
end

 clearvars('-except',main_preload_vars{:});
