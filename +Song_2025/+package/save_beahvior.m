%% reaction time , performance, iti move cross modality
clear all
Path = 'D:\Data process\project_cross_model\wf_data\';
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data';

behavior_aligned=cell(2,1);
behavior_each_mice=cell(2,1);

for curr_group=1:2
    switch curr_group
        case 1
            animals = {'DS007','DS010','AP019','AP021','DS011','AP022'};n1_name='visual position';n2_name='audio volume';
        case 2
            animals = {'DS000','DS004','DS014','DS015','DS016'};n1_name='audio volume';n2_name='visual position';
    end

    p_val=cell(length(animals),1);
    reaction_time=cell(length(animals),1);
    reward_time=cell(length(animals),1);
    reaction_time_null=cell(length(animals),1);

    workflow_name=cell(length(animals),1);
    performance=cell(length(animals),1);

    itimove=cell(length(animals),1);
    itimove_all=cell(length(animals),1);
velocity=cell(length(animals),1);
    for curr_animal =1:length(animals)
        animal=animals{curr_animal};
        raw_data_behavior=load([Path   'behavior\' animal '_behavior'  '.mat']);
        raw_data_task=load([Path   'task\' animal '_task'  '.mat']);

        if strcmp(animal,'AP019')
            temp_idx=  (1:15)';
        else
            [~, temp_idx] = ismember( raw_data_task.workflow_day,raw_data_behavior.workflow_day);
        end

        tem_p=nan(length(raw_data_behavior.workflow_day),2);
        idx_v=ismember(raw_data_behavior.workflow_name,'visual position');
        idx_a=ismember(raw_data_behavior.workflow_name,'audio volume');
        idx_m=ismember(raw_data_behavior.workflow_name,'mixed VA');

        tem_p(idx_v,1)= raw_data_behavior.rxn_l_mad_p(idx_v,1);
        tem_p(idx_a,1)=raw_data_behavior.rxn_l_mad_p(idx_a,1);
        tem_p(idx_m,:)= [raw_data_behavior.rxn_l_mad_p(idx_m,2)...
            raw_data_behavior.rxn_l_mad_p(idx_m,3)];
        % p_val{curr_animal}=tem_p(temp_idx)<0.01;

        tem_p_val=nan(length(raw_data_behavior.workflow_day),1);
        tem_p_val(idx_v,1)= raw_data_behavior.rxn_l_mad_p(idx_v,1)<0.01;
        tem_p_val(idx_a,1)=raw_data_behavior.rxn_l_mad_p(idx_a,1)<0.01;
        tem_p_val(idx_m,1)= raw_data_behavior.rxn_l_mad_p(idx_m,2)<0.01 &...
            raw_data_behavior.rxn_l_mad_p(idx_m,3)<0.01;

        p_val{curr_animal}=tem_p_val(temp_idx);


        temp_reaction_time=nan(length(raw_data_behavior.workflow_day),2);
        temp_reaction_time(idx_v)= raw_data_behavior.stim2move_med(idx_v,1);
        temp_reaction_time(idx_a)=raw_data_behavior.stim2move_med(idx_a,1);
        temp_reaction_time(idx_m,:)= [raw_data_behavior.stim2move_med(idx_m,2)...
            raw_data_behavior.stim2move_med(idx_m,3)];

        temp_reaction_time_null=nan(length(raw_data_behavior.workflow_day),2);
        temp_reaction_time_null(idx_v)= raw_data_behavior.stim2move_med_null(idx_v,1);
        temp_reaction_time_null(idx_a)=raw_data_behavior.stim2move_med_null(idx_a,1);
        temp_reaction_time_null(idx_m,:)= [raw_data_behavior.stim2move_med_null(idx_m,2)...
            raw_data_behavior.stim2move_med_null(idx_m,3)];


        reaction_time{curr_animal}=temp_reaction_time(temp_idx,:);
        reaction_time_null{curr_animal}=temp_reaction_time_null(temp_idx,:);


        temp_reaction_time_mad=nan(length(raw_data_behavior.workflow_day),2);
        temp_reaction_time_mad(idx_v)= raw_data_behavior.stim2lastmove_mad(idx_v,1);
        temp_reaction_time_mad(idx_a)=raw_data_behavior.stim2lastmove_mad(idx_a,1);
        temp_reaction_time_mad(idx_m,:)= [raw_data_behavior.stim2lastmove_mad(idx_m,2)...
            raw_data_behavior.stim2lastmove_mad(idx_m,3)];

        temp_reaction_time_mad_null=nan(length(raw_data_behavior.workflow_day),2);
        temp_reaction_time_mad_null(idx_v)= raw_data_behavior.stim2lastmove_mad_null(idx_v,1);
        temp_reaction_time_mad_null(idx_a)=raw_data_behavior.stim2lastmove_mad_null(idx_a,1);
        temp_reaction_time_mad_null(idx_m,:)= [raw_data_behavior.stim2lastmove_mad_null(idx_m,2)...
            raw_data_behavior.stim2lastmove_mad_null(idx_m,3)];

        performance{curr_animal}=(temp_reaction_time_mad_null(temp_idx,:)-temp_reaction_time_mad(temp_idx,:))./...
            (temp_reaction_time_mad(temp_idx,:)+temp_reaction_time_mad_null(temp_idx,:));
  
        

        temp_velocity=nan(length(raw_data_behavior.workflow_day),2);
        temp_velocity(idx_v)= cellfun(@(x)  -min(nanmean(x,1)) ,raw_data_behavior.frac_velocity_stimalign(idx_v,1),'UniformOutput',true);
        temp_velocity(idx_a)=cellfun(@(x) -min(nanmean(x,1)) ,raw_data_behavior.frac_velocity_stimalign(idx_a,1),'UniformOutput',true);
        if sum(idx_m)>0
        temp_velocity(idx_m,:)= [cellfun(@(x) -min( nanmean(x,1)) ,raw_data_behavior.frac_velocity_stimalign(idx_m,2),'UniformOutput',true)...
           cellfun(@(x) -min( nanmean(x,1)) ,raw_data_behavior.frac_velocity_stimalign(idx_m,3),'UniformOutput',true)];
        end

        velocity{curr_animal}=temp_velocity;

      

        reward_time{curr_animal}=cellfun(@mean ,raw_data_behavior.stim_on2off_times(temp_idx),'UniformOutput',true);
        workflow_name{curr_animal}=raw_data_behavior.workflow_name(temp_idx) ;

        itimove{curr_animal}=cellfun(@(x,y) length(x)/length(y),raw_data_behavior.iti_move(temp_idx),  raw_data_behavior.stim2move_times((temp_idx),1),'UniformOutput',true);
        itimove_all{curr_animal}=cellfun(@(x,y) length(x)/length(y),raw_data_behavior.all_iti_move(temp_idx),  raw_data_behavior.stim2move_times((temp_idx),1),'UniformOutput',true);


    end
    behavior_each_mice{curr_group}.velocity=velocity;
    behavior_each_mice{curr_group}.reaction_time=reaction_time;
    behavior_each_mice{curr_group}.performance=performance;
    behavior_each_mice{curr_group}.itimove=itimove;
    behavior_each_mice{curr_group}.learned=p_val;
       behavior_each_mice{curr_group}.workflow_name=workflow_name;

    names={'reaction_time','performance','itimove','velocity'};
    for curr_dat=1:4
        behavior_aligned{curr_group}.(names{curr_dat})=table;
        behavior_aligned{curr_group}.(names{curr_dat}).name=animals';

        % naive
        performance_pre0=cell(length(animals{curr_group}),1);
        performance_pre0=cellfun(@(perf,p,name)  perf(find( strcmp(n1_name,name) &p(:,1)==0),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        performance_pre0 = cellfun(@(x) reshape(x(1:max(end-2,0)),[],1),performance_pre0,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod1_naive=performance_pre0;
        
        % pre_learn
        performance_pre1=cell(length(animals{curr_group}),1);
        performance_pre1=cellfun(@(perf,p,name)  perf(find( strcmp(n1_name,name) &p(:,1)==0,2,'last'),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod1_pre_learn=performance_pre1;


        % post_learn
        performance_post1=cell(length(animals),1);
        performance_post1=cellfun(@(perf,p,name)  perf(find( strcmp(n1_name,name) &p(:,1)==1,5,'first'),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod1_post_learn=performance_post1;

        % mod2_pre_learn
        performance_pre2=cell(length(animals),1);
        performance_pre2=cellfun(@(perf,p,name)  perf(find( strcmp(n2_name,name) &p(:,1)==0,2,'first'),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod2_pre_learn=performance_pre2;

        % mod2_post_learn
        performance_post2=cell(length(animals),1);
        performance_post2=cellfun(@(perf,p,name)  perf(find( strcmp(n2_name,name) &p(:,1)==1,5,'first'),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod2_post_learn=performance_post2;


        % mod2
        performance_first6day=cell(length(animals{curr_group}),1);
        performance_first6day = cellfun(@(perf,p,name)  perf(find( strcmp(n2_name,name),6,'first'),1 )  ,...
            behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
        behavior_aligned{curr_group}.(names{curr_dat}).mod2=performance_first6day;

        % mix
        switch curr_dat
            case {1,2,4}
                performance_post3_v=cell(length(animals),1);
                performance_post3_v=cellfun(@(perf,p,name)  perf(find( strcmp('mixed VA',name) &p(:,1)==1,3,'first'),1 )  ,...
                    behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);
               
                performance_post3_a=cell(length(animals),1);

                performance_post3_a=cellfun(@(perf,p,name)  perf(find( strcmp('mixed VA',name) &p(:,1)==1,3,'first'),2 )  ,...
                    behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);

                behavior_aligned{curr_group}.(names{curr_dat}).mix_v=performance_post3_v;
                behavior_aligned{curr_group}.(names{curr_dat}).mix_a=performance_post3_a;
            case 3
                performance_post3=cell(length(animals),1);
                performance_post3=cellfun(@(perf,p,name)  perf(find( strcmp('mixed VA',name) &p(:,1)==1,3,'first'),1 )  ,...
                    behavior_each_mice{curr_group}.(names{curr_dat}),p_val,workflow_name ,'UniformOutput',false);

                behavior_aligned{curr_group}.(names{curr_dat}).mix=performance_post3;
        end
    end
end


% filled with nan matrix
stage_dim=[nan,1 2 5 1 5 6 3 3];
tmp = cellfun(@(group) structfun(@(wk)  arrayfun(@(a) arrayfun(@(s) ...
    vertcat( wk{a,s}{1}, ...
    repmat(NaN, ...
    max(stage_dim(s) - size(wk{a,s}{1},1)), 1) ), ...
    2:size(wk,2), 'UniformOutput', false), ...
    1:size(wk,1), 'UniformOutput', false), group,'UniformOutput',false),    behavior_aligned,'UniformOutput',false);


temp_data =cellfun(@(group) structfun(@(x)  vertcat(x{:}), group,'UniformOutput',false), tmp,'UniformOutput',false);   % 现在就是 6×10 的 cell，等价于上面循环版的结果
temp_data= cellfun(@(group) structfun(@(x)  subsasgn(x,substruct('()',{':',[1 4]}), ...
    cellfun(@(a) nanmean(a),x(:,[1 4]),'UniformOutput',false)), group,'UniformOutput',false), ...
    temp_data, 'UniformOutput', false);

temp_1= cellfun(@(group) structfun(@(x)  arrayfun(@(id) vertcat(x{id,:}),1:size(x,1),'UniformOutput',false),...
    group,'UniformOutput',false),temp_data,'UniformOutput',false);


behavior_across_day=  cellfun(@(group) ...
    structfun(@(fld) cat(2, fld{:}), group, 'UniformOutput', false), ...
    temp_1, 'UniformOutput', false);

save(fullfile(save_path,'behavior.mat' ),'behavior_each_mice','behavior_aligned','behavior_across_day','-v7.3');
