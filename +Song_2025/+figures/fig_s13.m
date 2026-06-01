%% EDF 13 audio balance to visual
main_preload_vars = who;

animals = {'DS029','DS030','DS031'};

temp_data_all=table;
for curr_animal=1:length(animals)
    preload_vars=who;

    animal=animals{curr_animal};
    data_all=load(fullfile(Path,'data\revision\single_data',[animal '_all_data.mat']));


    select_id_0=find((strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume_earphone_balance')|...
        strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume_earphone_balance')));
    select_id_0=select_id_0(cellfun(@(x)  x.rxn_l_p(1)<0.05 , data_all.behavior_task(select_id_0),'uni',true)==0);

    all_audio_task_pre=cellfun(@(x) x.stim_kernels,data_all.wf_task(select_id_0),'UniformOutput',false);



    select_id_1=find((strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume_earphone_balance')|...
        strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume_earphone_balance'))&...
        ~cellfun(@isempty ,data_all.wf_hml_passive_audio_earphone_balance_only),1,'first');
    % select_id_1=select_id_1(cellfun(@(x)  x.rxn_l_p(1)<0.05 , data_all.behavior_task(select_id_1),'uni',true)==0);

    all_hml_passive_pre=cellfun(@(x) cat(3,x.kernels_decoding{:}),data_all.wf_hml_passive_audio_earphone_balance_only(select_id_1),'UniformOutput',false);





    select_id_2=find((strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume_earphone_balance')|...
        strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume_earphone_balance'))&...
        ~cellfun(@isempty ,data_all.wf_hml_passive_audio_earphone_balance_only),2,'last');

    all_hml_passive_post=cellfun(@(x) cat(3,x.kernels_decoding{:}),data_all.wf_hml_passive_audio_earphone_balance_only(select_id_2),'UniformOutput',false);
    all_audio_task_post=cellfun(@(x) x.stim_kernels,data_all.wf_task(select_id_2),'UniformOutput',false);


    select_id_3=find((strcmp([data_all.task_name],'stim_wheel_right_stage1')|strcmp([data_all.task_name],'stim_wheel_right_stage2'))...
        &~cellfun(@isempty ,data_all.wf_lcr_passive),2,'last');

    all_visual_task=cellfun(@(x) x.stim_kernels,data_all.wf_task(select_id_3),'UniformOutput',false);

    all_lcr_passive=cellfun(@(x) cat(3,x.kernels_decoding{:}),data_all.wf_lcr_passive(select_id_3),'UniformOutput',false);
    % all_task=cellfun(@(x) cat(3,x.kernels_decoding{:}),data_all.wf_lcr_passive(select_id),'UniformOutput',false);

    temp_data_all.all_hml_passive_pre{curr_animal}=all_hml_passive_pre;
    temp_data_all.all_hml_passive_post{curr_animal}=all_hml_passive_post;

    temp_data_all.all_audio_task_pre{curr_animal}=all_audio_task_pre;
    temp_data_all.all_audio_task_post{curr_animal}=all_audio_task_post;

    temp_data_all.all_lcr_passive{curr_animal}=all_lcr_passive;
    temp_data_all.all_visual_task{curr_animal}=all_visual_task;
    clearvars('-except',preload_vars{:});
end

temp_hml_passive_pre= feval(@(a) cat(4,a{:}) ,cellfun(@(x) nanmean(cat(4,x{:}),4),temp_data_all.all_hml_passive_pre,'UniformOutput',false  ));
all_passive_image_pre= plab.wf.svd2px(U_master(:,:,1:size(cat(3,temp_hml_passive_pre),1)),temp_hml_passive_pre);
all_audio_passive_image_max_pre=permute(max(nanmean(all_passive_image_pre(:,:,kernels_period,:,:),5),[],3),[1,2,4,3]);
temp_wf_passive_plot_tace= ds.make_each_roi(all_passive_image_pre, length(t_kernels),roi1);


temp_hml_passive_post= feval(@(a) cat(4,a{:}) ,cellfun(@(x) nanmean(cat(4,x{:}),4),temp_data_all.all_hml_passive_post,'UniformOutput',false  ));
all_passive_image_post= plab.wf.svd2px(U_master(:,:,1:size(cat(3,temp_hml_passive_post),1)),temp_hml_passive_post);
all_audio_passive_image_max_post=permute(max(nanmean(all_passive_image_post(:,:,kernels_period,:,:),5),[],3),[1,2,4,3]);
temp_wf_passive_plot_tace= ds.make_each_roi(all_passive_image_pre, length(t_kernels),roi1);


temp_lcr_passive= feval(@(a) cat(4,a{:}) ,cellfun(@(x) nanmean(cat(4,x{:}),4),temp_data_all.all_lcr_passive,'UniformOutput',false  ));
all_visual_passive_image= plab.wf.svd2px(U_master(:,:,1:size(cat(3,temp_lcr_passive),1)),temp_lcr_passive);
all_visual_passive_image_max=permute(max(nanmean(all_visual_passive_image(:,:,kernels_period,:),5),[],3),[1,2,4,3]);
temp_visual_passive_plot_tace= ds.make_each_roi(all_visual_passive_image, length(t_kernels),roi1);

temp_visual=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_visual_task{:})));
all_visual_task_image= plab.wf.svd2px(U_master(:,:,1:size(temp_visual,1)),temp_visual);
all_visual_task_image_max=max(nanmean(all_visual_task_image(:,:,kernels_period,:),4),[],3);

temp_audio_pre=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_audio_task_pre{:})));
all_audio_task_image_pre= plab.wf.svd2px(U_master(:,:,1:size(temp_audio_pre,1)),temp_audio_pre);
all_audio_task_image_max_pre=max(nanmean(all_audio_task_image_pre(:,:,kernels_period,:),4),[],3);


temp_audio_post=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_audio_task_post{:})));
all_audio_task_image_post= plab.wf.svd2px(U_master(:,:,1:size(temp_audio_post,1)),temp_audio_post);
all_audio_task_image_max_post=max(nanmean(all_audio_task_image_post(:,:,kernels_period,:),4),[],3);


figure('Position',[50 50 600 300]);
tiledlayout(2,3,'TileIndexing','columnmajor')
colors=[[0 0 1];[0 0 0];[ 1 0 0]];
used_passive=[3,3,3];
names={'A-pre-learning','A-well-trained','V-well-trained'}
color_n={'R','R','B'}
for curr_stage=1:3
    switch curr_stage
        case 1
            tem_image=all_audio_passive_image_max_pre;
            tem_trace_mean=permute(nanmean(temp_wf_passive_plot_tace(3,:,:,:),4),[2,3,1]);
            tem_trace_error=permute(nanstd(temp_wf_passive_plot_tace(3,:,:,:),0,4)./sqrt(size(temp_wf_passive_plot_tace,4)),[2,3,1]);
            temp_task=all_audio_task_image_max_pre;
        case 2

            tem_image=all_audio_passive_image_max_post;
            tem_trace_mean=permute(nanmean(temp_wf_passive_plot_tace(3,:,:,:),4),[2,3,1]);
            tem_trace_error=permute(nanstd(temp_wf_passive_plot_tace(3,:,:,:),0,4)./sqrt(size(temp_wf_passive_plot_tace,4)),[2,3,1]);
            temp_task=all_audio_task_image_max_post;
        case 3
            tem_image=all_visual_passive_image_max;
            tem_trace_mean=permute(nanmean(temp_visual_passive_plot_tace(1,:,:,:),4),[2,3,1]);
            tem_trace_error=permute(nanstd(temp_visual_passive_plot_tace(1,:,:,:),0,4)./sqrt(size(temp_visual_passive_plot_tace,4)),[2,3,1]);

            temp_task=all_visual_task_image_max;


    end
    a1=nexttile
    imagesc(temp_task);
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim(0.0004 .* [ 0, 1]);
    colormap(a1, ap.colormap(['W'  color_n{curr_stage}] ));
    title(names{curr_stage},'FontWeight','normal')

    for curr_passive=used_passive(curr_stage)
        a2=nexttile
        imagesc(tem_image(:,:,curr_passive));
        axis image off;
        ap.wf_draw('ccf', [0.5 0.5 0.5]);
        clim(0.0004 .* [ 0, 1]);
        colormap(a2, ap.colormap(['W'  color_n{curr_stage}] ));
    end
    colorbar

end
colorbar


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_D.eps'), ...
    'ContentType','vector');

 clearvars('-except',main_preload_vars{:});
