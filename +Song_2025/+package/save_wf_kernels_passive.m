clear all
clc
Path = 'D:\Data process\wf_data\';
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data';

all_workflow={'lcr_passive','hml_passive_audio'};
Color={'B','R'};

n1_name='';n2_name='';
animals={};
groups={'V_A','A_V'};

data_all_days=cell(2,1);
data_type={'raw','kernels'};


surround_samplerate = 35;
surround_window_passive = [-0.5,1];
surround_window_task = [-0.2,1];
t_kernels=1/surround_samplerate*[-10:30];
passive_boundary=0.2;
period_kernels=find(t_kernels>0&t_kernels<passive_boundary);



wf_passive_kernel_aligned=struct;
wf_passive_kernel_each_mice=cell(2,1);

for  workflow_idx=1:2
    workflow=all_workflow{workflow_idx};
    wf_passive_kernel_aligned.(workflow)=cell(2,1);
    wf_passive_kernel_each_mice{workflow_idx}=cell(2,1);

    for curr_group=1:2
        main_preload_vars = who;
        if curr_group==1
            animals{curr_group} = {'DS007','DS010','AP019','AP021','DS011','AP022'};n1_name='visual position';n2_name='audio volume';
        elseif curr_group==2
            animals{curr_group} = {'DS000','DS004','DS014','DS015','DS016'};n1_name='audio volume';n2_name='visual position';

        end

        wf_passive_kernel_aligned.(workflow){curr_group}=table;
        wf_passive_kernel_aligned.(workflow){curr_group}.name= animals{curr_group}';

        % used_id=1:3;
        all_data_video=cell(length(animals{curr_group}),1);
        all_data_workflow_name=cell(length(animals{curr_group}),1);
        all_data_learned_day=cell(length(animals{curr_group}),1);
        matches=cell(length(animals{curr_group}),1);

        for curr_animal=1:length(animals{curr_group})
            preload_vars = who;
            animal=animals{curr_group}{curr_animal};
            raw_data_passive=load([Path  workflow '\' animal '_' workflow '.mat']);
            raw_data_behavior=load([Path   'behavior\' animal '_behavior'  '.mat']);


            idx_all=cellfun(@(x) ~isempty(x),raw_data_passive.wf_px_kernels);
            matches{curr_animal}=unique(raw_data_passive.workflow_type_name_merge(idx_all)  ,'stable');

            all_data_video{curr_animal}=raw_data_passive.wf_px_kernels(idx_all)';

            all_data_workflow_name{curr_animal}=raw_data_passive.workflow_type_name_merge(idx_all);


            [~, temp_idx] = ismember( raw_data_passive.workflow_day,raw_data_behavior.workflow_day);

            temp_p=nan(length(raw_data_behavior.workflow_day),2);
            idx_v=ismember(raw_data_behavior.workflow_name,'visual position');
            idx_a=ismember(raw_data_behavior.workflow_name,'audio volume');
            idx_m=ismember(raw_data_behavior.workflow_name,'mixed VA');

            temp_p(idx_v,1)= raw_data_behavior.rxn_l_mad_p(idx_v,1);
            temp_p(idx_a,1)=raw_data_behavior.rxn_l_mad_p(idx_a,1);
            temp_p(idx_m,:)= [raw_data_behavior.rxn_l_mad_p(idx_m,1)...
                raw_data_behavior.rxn_l_mad_p(idx_m,2)];
            temp_p=temp_p(temp_idx(temp_idx>0),:);


            if workflow_idx==1
                temp_p1= [nan(sum(temp_idx==0 ),1) ;  temp_p(isnan(temp_p(:,2)),1) ;   temp_p(~isnan(temp_p(:,2)),1)];
            else
                temp_p1= [nan(sum(temp_idx==0 ),1) ;  temp_p(isnan(temp_p(:,2)),1) ;   temp_p(~isnan(temp_p(:,2)),2)];
            end

            all_data_learned_day{curr_animal}=temp_p1(idx_all)<0.01;

            clearvars('-except',preload_vars{:});

        end

        wf_passive_kernel_each_mice{curr_group}{workflow_idx}=all_data_video;
        % habituation
        habituation_idx=cellfun(@(x) any(strcmp('naive',x)),matches,'UniformOutput',true );
        habituation_data=cell(length(animals{curr_group}),1);
        habituation_data(habituation_idx) =  cellfun(@(x, y, z)...
            x(find(strcmp(y, z(find(cellfun(@(idx) strcmp('naive', idx), z, 'UniformOutput', true)))), 3, 'first'), :), ...
            all_data_video(habituation_idx), all_data_workflow_name(habituation_idx), matches(habituation_idx), ...
            'UniformOutput', false);

        wf_passive_kernel_aligned.(workflow){curr_group}.habituation=habituation_data;

      
        % mod1_naive

        pre_learn_data0=cell(length(animals{curr_group}),1);
        pre_learn_data0 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l==0 ))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        pre_learn_data0 = cellfun(@(x) reshape(x(1:max(end-2,0)),[],1),pre_learn_data0,'UniformOutput',false);

        wf_passive_kernel_aligned.(workflow){curr_group}.mod1_naive=pre_learn_data0;


        % mod1_pre_learn
        pre_learn_data1=cell(length(animals{curr_group}),1);
        pre_learn_data1 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l==0 ,2,'last'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        wf_passive_kernel_aligned.(workflow){curr_group}.mod1_pre_learn=pre_learn_data1;


        % mod1_post_learn
        post_learn1_data1=cell(length(animals{curr_group}),1);
        post_learn1_data1 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l==1 ,5,'first'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);

        wf_passive_kernel_aligned.(workflow){curr_group}.mod1_post_learn=post_learn1_data1;


        % 
        % post_learn2_data1=cell(length(animals{curr_group}),1);
        % post_learn2_data1 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l==1 ,5,'first'))...
        %     ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        % post_learn2_data1 = cellfun(@(x) x(4:end),post_learn2_data1,'UniformOutput',false);
        % 
        % wf_passive_kernel_aligned.(workflow){curr_group}.mod1_well_trained=post_learn2_data1;



      

        %mod2_pre_learn
        pre_learn_data2=cell(length(animals{curr_group}),1);
        pre_learn_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l==0 ,2,'first'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);

        wf_passive_kernel_aligned.(workflow){curr_group}.mod2_pre_learn=pre_learn_data2;


        %mod2_post_learn
        post_learn1_data2=cell(length(animals{curr_group}),1);
        post_learn1_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l==1 ,5,'first'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        wf_passive_kernel_aligned.(workflow){curr_group}.mod2_post_learn=post_learn1_data2;


        % post_learn2_data2=cell(length(animals{curr_group}),1);
        % post_learn2_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true))))& l==1 ,5,'first'))...
        %     ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        % post_learn2_data2 = cellfun(@(x) x(4:end),post_learn2_data2,'UniformOutput',false);
        % wf_passive_kernel_aligned.(workflow){curr_group}.mod2_well_trained=post_learn2_data2;


        %mod2
        first_6day_data2=cell(length(animals{curr_group}),1);
        first_6day_data2 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n2_name, idx),z,'UniformOutput',true)))) ,6,'first'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        wf_passive_kernel_aligned.(workflow){curr_group}.mod2=first_6day_data2;


        n3_name='mixed VA';
        mixed_idx=cellfun(@(x) any(strcmp(n3_name, x)),matches ,'UniformOutput',true);
        data3=cell(length(animals{curr_group}),1);
        data3(mixed_idx) = cellfun(@(x,y,z) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n3_name, idx),z,'UniformOutput',true)))),3,'last'),:),...
            all_data_video(mixed_idx) ,all_data_workflow_name(mixed_idx) ,matches(mixed_idx) ,'UniformOutput',false);

        wf_passive_kernel_aligned.(workflow){curr_group}.mix=data3;


        clearvars('-except',main_preload_vars{:});

    end
end



% filled with nan matrix
stage_dim=[nan,3 1 2 5 1 5 6 3 ];
tmp =structfun(@(A) cellfun(@(wk)   arrayfun(@(a) arrayfun(@(s) ...
    vertcat( wk{a,s}{1}, ...
    repmat({NaN(200, numel(t_kernels),3)}, ...
    max(stage_dim(s) - size(wk{a,s}{1},1), 0), 1) ), ...
    2:9, 'UniformOutput', false), ...
    1:size(wk,1), 'UniformOutput', false), A,'UniformOutput',false),wf_passive_kernel_aligned,'UniformOutput',false);

% wk=wf_passive_kernel_aligned.lcr_passive{1}

temp_data = structfun(@(a) cellfun(@(x) vertcat(x{:}),a,'UniformOutput',false),tmp,'UniformOutput',false);   % 现在就是 6×10 的 cell，等价于上面循环版的结果
temp_data= structfun(@(A) cellfun(@(x) subsasgn(x,substruct('()',{':',[2 5]}), ...
    cellfun(@(a) {nanmean(cat(4,a{:}),4)},x(:,[2 5]),'UniformOutput',false)),A,'UniformOutput',false), ...
    temp_data, 'UniformOutput', false);

min_dim=200;
temp_data=structfun(@(a) cellfun(@(x) cellfun(@(n) cellfun(@(a)   a(1:min_dim,:,:),n,'UniformOutput',false),x,'UniformOutput',false ),a,'UniformOutput',false),temp_data, 'UniformOutput', false );
temp_1= structfun(@(a)  cellfun(@(x)  arrayfun(@(id) vertcat(x{id,:}),1:size(x,1),'UniformOutput',false   ),a,'UniformOutput',false),temp_data,'UniformOutput',false);
temp_2= structfun(@(a) cellfun(@(x)    feval(@(C) cat(5, C{:}),cellfun(@(a)   cat(4,a{:}),x,'UniformOutput',false)),a,'UniformOutput',false),temp_1,'UniformOutput',false);
wf_passive_kernels_across_day ={temp_2.lcr_passive;temp_2.hml_passive_audio};

save(fullfile(save_path, 'wf_passive_kernels.mat' ),'wf_passive_kernel_each_mice','wf_passive_kernel_aligned','wf_passive_kernels_across_day','-v7.3')
