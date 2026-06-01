%% Generate figures for Song et al 2025
clear all
clc
Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';
U_master = plab.wf.load_master_U;
load(fullfile(Path,'data\General_information\roi.mat'))


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

for curr_figure =1:6
    fprintf('Starting drawing Figure %d...\n', curr_figure);
    feval(sprintf('Song_2025.figures.figure_%d', curr_figure));
    fprintf('Finished Figure %d.\n', curr_figure);
end

for curr_figs = 1:18
    fprintf('Starting drawing Figure S%d...\n', curr_figs);
    feval(sprintf('Song_2025.figures.fig_s%d', curr_figs));
    fprintf('Finished Figure S%d.\n', curr_figs);
end













