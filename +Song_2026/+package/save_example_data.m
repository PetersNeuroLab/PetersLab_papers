%% fig s1 s2
clear all
Path = 'D:\Data process\wf_data\';
save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data';


%% fig s1a  example trace 

main_preload_vars = who;
animal ='DS019';
rec_day='2025-01-23';
rec = plab.find_recordings(animal,rec_day,'*wheel*');
rec_time = rec.recording{end};
load_parts = struct;
load_parts.behavior = true;
load_parts.widefield = false;
ap.load_recording;
new_vars = setdiff(who, main_preload_vars);

save(fullfile(save_path, 'example_trace.mat'),new_vars{:});

clearvars('-except',main_preload_vars{:});


%% fig s1b

main_preload_vars = who;
animals = {'DS010','DS015'};
behavior=struct;
for curr_animal=1:2
animal=animals{curr_animal};
behavior.(animal)=load([Path   'behavior\' animal '_behavior'  '.mat']);
end
save(fullfile(save_path, 'example_behaviors.mat'),'behavior','-v7.3');
clearvars('-except',main_preload_vars{:});


%% fig s2d passive_average
main_preload_vars = who;
all_workflow={'lcr_passive','hml_passive_audio'};

wf_passive_average=struct;
for  workflow_idx=1:2
    workflow=all_workflow{workflow_idx};
    wf_passive_average.(workflow)=cell(2,1);

    for curr_group=workflow_idx
        preload_vars = who;
        if curr_group==1
            animals{curr_group} = {'DS007','DS010','AP019','AP021','DS011','AP022'};n1_name='visual position';n2_name='audio volume';
        elseif curr_group==2
            animals{curr_group} = {'DS000','DS004','DS014','DS015','DS016'};n1_name='audio volume';n2_name='visual position';

        end

        % wf_passive_kernel_aligned.(workflow){curr_group}=table;
        % wf_passive_kernel_aligned.(workflow){curr_group}.name= animals{curr_group}';

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

            all_data_video{curr_animal}=raw_data_passive.wf_px(idx_all)';

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
      
        % mod1_post_learn
        post_learn1_data1=cell(length(animals{curr_group}),1);
        post_learn1_data1 = cellfun(@(x,y,z,l) x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l==1 ,5,'first'))...
            ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
        
        post_learn1_data1=cellfun(@(x)  nanmean(cat(4,x{min(4,length(x)):end}),4) ,post_learn1_data1,'UniformOutput',false );

        wf_passive_average.(workflow)=post_learn1_data1;

        clearvars('-except',preload_vars{:});

    end
end

save(fullfile(save_path, 'wf_passive_average.mat' ),'wf_passive_average','-v7.3')
clearvars('-except',main_preload_vars{:});

%% fig s2a task_average
main_preload_vars = who;
workflow='task';
wf_task_average_aligned=cell(2,1);
wf_task_encoding_aligned=cell(2,1);

for curr_group=1:2
    main_preload_vars = who;
    if curr_group==1
        animals{curr_group} = {'DS007','DS010','AP019','AP021','DS011','AP022'};n1_name='visual position';n2_name='audio volume';
    elseif curr_group==2
        animals{curr_group} = {'DS000','DS004','DS014','DS015','DS016'};n1_name='audio volume';n2_name='visual position';

    end
    all_data_video=cell(length(animals{curr_group}),1);
        all_data_encoding=cell(length(animals{curr_group}),1);

    all_data_workflow_name=cell(length(animals{curr_group}),1);
    all_data_learned_day=cell(length(animals{curr_group}),1);
    matches=cell(length(animals{curr_group}),1);

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

        idx=cellfun(@(x) ~isempty(x),raw_data_task.wf_px_task);
        
        image_all_average(idx)=cellfun(@(x)  x(:,:,1),raw_data_task.wf_px_task(idx),'UniformOutput',false);
        image_all_encoding(idx)=cellfun(@(x)  permute(x(1,:,:),[3,2,1]) ,raw_data_task.wf_px_task_kernels_encode(idx),'UniformOutput',false);

        matches{curr_animal}=unique(raw_data_task.workflow_type_name_merge(idx)  ,'stable');

        all_data_video{curr_animal}=image_all_average(idx)';
        all_data_encoding{curr_animal}=image_all_encoding(idx)';

        all_data_workflow_name{curr_animal}=raw_data_task.workflow_type_name_merge(idx);
        all_data_learned_day{curr_animal}=raw_data_task.learned_day(idx,:);
        clearvars('-except',preload_vars{:});

    end


    % mod1_post_learn
    post_learn1_data1=cell(length(animals{curr_group}),1);
    post_learn1_data1 = cellfun(@(x,y,z,l)...
        x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
        ,all_data_video,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    post_learn1_data1=cellfun(@(x)   nanmean(cat(4,x{min(length(x),4):end}),4)  ,post_learn1_data1,'UniformOutput',false);
    wf_task_average_aligned{curr_group}=post_learn1_data1;

      % mod1_post_learn
    post_learn1_data1_encoding=cell(length(animals{curr_group}),1);
    post_learn1_data1_encoding = cellfun(@(x,y,z,l)...
        x(find(strcmp(y,z(find(cellfun(@(idx) strcmp(n1_name, idx),z,'UniformOutput',true))))& l(:,1)==1 ,5,'first'))...
        ,all_data_encoding,all_data_workflow_name,matches,all_data_learned_day,'UniformOutput',false);
    post_learn1_data1_encoding=cellfun(@(x)   nanmean(cat(4,x{min(length(x),4):end}),4)  ,post_learn1_data1_encoding,'UniformOutput',false);
    
    
    wf_task_encoding_aligned{curr_group}=post_learn1_data1_encoding;


    clearvars('-except',main_preload_vars{:});
end

save(fullfile(save_path,'wf_task_average_encoding.mat' ),'wf_task_average_aligned','wf_task_encoding_aligned','-v7.3')
clearvars('-except',main_preload_vars{:});


