clear all
clc
Path = 'D:\Data process\wf_data\';
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data';

surround_samplerate = 35;
surround_window_task = [-0.2,1];
t_kernels=1/surround_samplerate*[-10:30];
task_boundary1=0;
task_boundary2=0.2;
period_kernels=find(t_kernels>task_boundary1&t_kernels<task_boundary2);

n1_name='';n2_name='';
use_period=[];
use_t=[];
animals={};
groups={'V_A','A_V'};
workflow='task';

wf_task_kernel_aligned=cell(2,1);
wf_task_kernel_each_mice=cell(2,1);
for curr_group=1:2;
    wf_task_kernel_aligned{curr_group}=table;
    main_preload_vars = who;
    if curr_group==1
        animals{curr_group} = {'DS007','DS010','AP019','AP021','DS011','AP022'};n1_name='visual position';n2_name='audio volume';
    elseif curr_group==2
        animals{curr_group} = {'DS000','DS004','DS014','DS015','DS016'};n1_name='audio volume';n2_name='visual position';

    end
    wf_task_kernel_aligned{curr_group}.name=animals{curr_group}';
    all_data_video=cell(length(animals{curr_group}),1);
    all_data_workflow_name=cell(length(animals{curr_group}),1);
    all_data_learned_day=cell(length(animals{curr_group}),1);
    matches=cell(length(animals{curr_group}),1);
    % all_data_sim2move=cell(length(animals{curr_group}),1);
    use_t=[];
    use_period=[];

    for curr_animal=1:length(animals{curr_group})
        preload_vars = who;

        animal=animals{curr_group}{curr_animal};
        raw_data_task=load([Path  workflow '\' animal '_' workflow '.mat']);
        raw_data_behavior=load([Path   'behavior\' animal '_behavior'  '.mat']);
        [~, temp_idx] = ismember( raw_data_task.workflow_day,raw_data_behavior.workflow_day);
        temp_p=nan(length(raw_data_behavior.workflow_day),2);
        idx_v=ismember(raw_data_behavior.workflow_name,'visual position');
        idx_a=ismember(raw_data_behavior.workflow_name,{'audio volume','audio freque'});
        idx_m=ismember(raw_data_behavior.workflow_name,'mixed VA');
        temp_p(idx_v,1)= raw_data_behavior.rxn_l_mad_p(idx_v,1);
        temp_p(idx_a,1)=raw_data_behavior.rxn_l_mad_p(idx_a,1);
        temp_p(idx_m,:)= [raw_data_behavior.rxn_l_mad_p(idx_m,1)...
            raw_data_behavior.rxn_l_mad_p(idx_m,2)];
        temp_p=temp_p(temp_idx,:);
        raw_data_task.learned_day=temp_p<0.01;
      
        idx=cellfun(@(x) ~isempty(x),raw_data_task.wf_px_task_kernels);
        image_all(idx)=cellfun(@(x)  x{1},raw_data_task.wf_px_task_kernels(idx),'UniformOutput',false);
        use_period=period_kernels;
        use_t=t_kernels;

        matches{curr_animal}=unique(raw_data_task.workflow_type_name_merge(idx)  ,'stable');

        all_data_video{curr_animal}=image_all(idx)';
        all_data_workflow_name{curr_animal}=raw_data_task.workflow_type_name_merge(idx);
        all_data_learned_day{curr_animal}=raw_data_task.learned_day(idx,:);
        clearvars('-except',preload_vars{:});

    end

    wf_task_kernel_each_mice{curr_group}=all_data_video;

    % mod1_naive
    pre_learn_data0=cell(length(animals{curr_group}),1);
    pre_learn_data0 = cellfun(@(x,y,z,l) ...
        x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==0 ))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    pre_learn_data0 = cellfun(@(x) reshape(x(1:max(end-2,0)),[],1),pre_learn_data0,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod1_naive=pre_learn_data0;

    % mod1_pre_learn

    pre_learn_data1=cell(length(animals{curr_group}),1);
    pre_learn_data1 = cellfun(@(x,y,z,l) ...
        x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==0 ,2,'last'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod1_pre_learn=pre_learn_data1;


    % mod1_post_learn
    post_learn1_data1=cell(length(animals{curr_group}),1);
    post_learn1_data1 = cellfun(@(x,y,z,l)...
        x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod1_post_learn=post_learn1_data1;

    % % mod1_well_trained
    % post_learn2_data1=cell(length(animals{curr_group}),1);
    % post_learn2_data1 = cellfun(@(x,y,z,l)...
    %     x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
    %     ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    % post_learn2_data1 = cellfun(@(x) x(4:end),post_learn2_data1,'UniformOutput',false);
    % wf_task_kernel_aligned{curr_group}.mod1_well_trained=post_learn2_data1;


    % mod2_pre_learn
    pre_learn_data2=cell(length(animals{curr_group}),1);
    pre_learn_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l(:,1)==0 ,2,'first'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod2_pre_learn=pre_learn_data2;

    % mod2_post_learn
    post_learn1_data2=cell(length(animals{curr_group}),1);
    post_learn1_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod2_post_learn=post_learn1_data2;

    % 
    % % mod2_well_trained
    % post_learn2_data2=cell(length(animals{curr_group}),1);
    % post_learn2_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
    %     ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    % post_learn2_data2 = cellfun(@(x) x(4:end),post_learn2_data2,'UniformOutput',false);
    % wf_task_kernel_aligned{curr_group}.mod2_well_trained=post_learn2_data2;

    % mod2_first_6_days
    first6day_data2=cell(length(animals{curr_group}),1);
    first6day_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true)))) ,6,'first'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mod2=first6day_data2;


    % mix_v& mix_a
    n3_name='mixed VA';
    mixed_idx=cellfun(@(x) any(strcmp(n3_name, x)),matches ,'UniformOutput',true);
    data3=cell(length(animals{curr_group}),1);
    data3(mixed_idx) = cellfun(@(x,y,z) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n3_name, idx),z,'UniformOutput',true)))),3,'last')),...
        all_data_video(mixed_idx) ,all_data_workflow_name(mixed_idx) ,matches(mixed_idx) ,'UniformOutput',false);


    wf_task_kernel_aligned{curr_group}.mix_v(mixed_idx)=cellfun(@(x) cellfun(@(n)  n(:,:,1)   ,x,'UniformOutput',false) ,data3(mixed_idx) ,'UniformOutput',false);
    wf_task_kernel_aligned{curr_group}.mix_a(mixed_idx)=cellfun(@(x) cellfun(@(n)  n(:,:,2)   ,x,'UniformOutput',false) ,data3(mixed_idx) ,'UniformOutput',false);

    clearvars('-except',main_preload_vars{:});
end


% filled with nan matrix
stage_dim=[nan,1 2 5 1 5 6 3 3];
tmp = cellfun(@(wk)   arrayfun(@(a) arrayfun(@(s) ...
    vertcat( wk{a,s}{1}, ...
    repmat({NaN(200, numel(use_t))}, ...
    max(stage_dim(s) - size(wk{a,s}{1},1)), 1) ), ...
    2:9, 'UniformOutput', false), ...
    1:size(wk,1), 'UniformOutput', false), wf_task_kernel_aligned,'UniformOutput',false);

temp_data =cellfun(@(x) vertcat(x{:}),tmp,'UniformOutput',false);   % 现在就是 6×10 的 cell，等价于上面循环版的结果
temp_data= cellfun(@(x) subsasgn(x,substruct('()',{':',[1 4]}), ...
    cellfun(@(a) {nanmean(cat(3,a{:}),3)},x(:,[1 4]),'UniformOutput',false)), ...
    temp_data, 'UniformOutput', false);

min_dim=170;
temp_data=cellfun(@(x) cellfun(@(n) cellfun(@(a)   a(1:min_dim,:),n,'UniformOutput',false),x,'UniformOutput',false ),temp_data, 'UniformOutput', false );
temp_1= cellfun(@(x)  arrayfun(@(id) vertcat(x{id,:}),1:size(x,1),'UniformOutput',false   ),temp_data,'UniformOutput',false);
wf_task_kernels_across_day= cellfun(@(x)    feval(@(C) cat(4, C{:}),cellfun(@(a)   cat(3,a{:}),x,'UniformOutput',false))   ,temp_1,'UniformOutput',false);

save(fullfile(save_path,'wf_task_kernels.mat' ),'wf_task_kernel_each_mice','wf_task_kernel_aligned','wf_task_kernels_across_day','-v7.3')
