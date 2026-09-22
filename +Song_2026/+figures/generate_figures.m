%% Set general paths and parameters
clear all
clc

% Set paths
Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2026';
U_master = plab.wf.load_master_U;
load(fullfile(Path,'data\General_information\roi.mat'))

% Set parameters
surround_samplerate = 35;
surround_window_task = [-0.2,1];
task_boundary1=0;
task_boundary2=0.2;

t_kernels=1/surround_samplerate*[-10:30];
kernels_period=find(t_kernels>task_boundary1&t_kernels<task_boundary2);

t_task = surround_window_task(1):1/surround_samplerate:surround_window_task(2);
period_task=find(t_task>0&t_task<0.2);

surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

surround_window = [-0.5,1];
mousecam_framerate = 30;
face_time = surround_window(1):1/mousecam_framerate:surround_window(2);


%% Draw figures

% Get all figure functions
fig_code_dir = dir(fullfile(fileparts(which('Song_2026.figures.generate_figures')),'*.m'));
fig_fcns = string(setdiff(erase({fig_code_dir.name},'.m'),'generate_figures'));

preload_vars = who;
for curr_fig_fcn = fig_fcns
    % Draw figure
    fprintf('Starting drawing Figure %s...\n', curr_fig_fcn);
    Song_2026.figures.(curr_fig_fcn);

    % Clear fig-related variables
    fprintf('Finished drawing Figure %s...\n', curr_fig_fcn);
    clearvars('-except',preload_vars{:});
end













