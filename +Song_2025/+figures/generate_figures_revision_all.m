%% EDF A    reaction time matched kernels


Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';
load(fullfile(Path,'data','revision','wf_task_decoding_concatnate.mat'));
load(fullfile(Path,'data\General_information\roi.mat'))

U_master = plab.wf.load_master_U;


local_data_path= 'D:\Data process\project_cross_model\wf_data\';
animals =     { 'DS007','DS010','AP019','AP021','DS011','AP022'...
    'DS000','DS004','DS014','DS015','DS016'};
surround_samplerate = 35;
t_kernels=1/surround_samplerate*[-10:30];

surround_window = [-0.2,1];
t_task = surround_window(1):1/surround_samplerate:surround_window(2);

surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);


% different brain regions roi
figure('Position',[50 50 200 200]);
merged_max=roi1(1).data.mask  +roi1(3).data.mask*2+...
    roi1(7).data.mask*3+roi1(9).data.mask*4+roi1(14).data.mask*5;
imagesc(merged_max)

axis image off;
cmap = [1 1 1
    0.40 0.60 0.85   % 蓝
    0.45 0.75 0.45   % 绿
    0.95 0.75 0.35   % 橙黄
    0.90 0.45 0.45   % 红
    0.65 0.50 0.85   % 紫
    ];

colormap(cmap);
h=ap.wf_draw('ccf', [0.5 0.5 0.5]);
% ap.wf_draw('grid')
% ap.wf_draw('area')
exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_A1.eps'), ...
    'ContentType','vector');


colors = [
    0.23 0.30 0.75
    0.40 0.50 0.85
    0.60 0.70 0.90
    0.80 0.85 0.95
    0.90 0.93 0.98
    ];

groups={'VA','AV'};
passive_workflows={'lcr_passive','hml_passive_audio'};
passive_id={3,2};
task_name={'visual position','audio volume'}
used_roi= [7 10];
cmap = [
    0.3 0.3 0.8;
    0.8 0.3 0.3;
    0.3 0.8 0.3;
    0.75 0.7 0.3;
    ];% 或者 lines(5), turbo(5)
% Colors = mat2cell(cmap, ones(size(temp_vel,1),1), 3);
Colors = mat2cell(cmap, ones(4,1), 3);
%
figure('Position',[50 50 800 600])
mainLayout = tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'none');
for curr_group=1:2
    animals_in_groups=wf_stim_kernels_concat.name(ismember(wf_stim_kernels_concat.group,groups{curr_group}));
    for curr_animal=1:length(animals_in_groups)
        animal=animals_in_groups{curr_animal};
        temp_passive_data=load(fullfile (local_data_path,passive_workflows{curr_group},[animal '_' passive_workflows{curr_group},'.mat']));
        temp_passive_wf= nanmean(cat(4,temp_passive_data.wf_px_kernels{find(  ismember (temp_passive_data.workflow_type_name_merge,task_name{curr_group}),5,'last')}),4);
        wf_stim_kernels_concat.wf_passive_kernels(find(contains(wf_stim_kernels_concat.name,animal)))={temp_passive_wf};
    end

    temp_wf_passive=wf_stim_kernels_concat.wf_passive_kernels(ismember(wf_stim_kernels_concat.group,groups{curr_group}));
    tem_image_passive=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),temp_wf_passive,'UniformOutput',false);
    tem_trace_passive=cellfun(@(x)   ds.make_each_roi(x,t_kernels,roi1),tem_image_passive,'UniformOutput',false);
    tem_trace_passive_mean=nanmean(cat(4,tem_trace_passive{:}),4);
    tem_trace_passive_error=std(cat(4,tem_trace_passive{:}),0,4,'omitmissing')./sqrt(length(tem_trace_passive));
    temp_wf=wf_stim_kernels_concat.wf_kernels(ismember(wf_stim_kernels_concat.group,groups{curr_group}));
    tem_image= feval(@(c) cat(2,c{:}),  cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
        ,temp_wf,'UniformOutput',false));
    temp_stim2move=feval(@(c) cat(2,c{:}), wf_stim_kernels_concat.stim_to_move(ismember(wf_stim_kernels_concat.group,groups{curr_group})));
    temp_stim2move_mean= nanmean(cellfun(@mean,temp_stim2move,'UniformOutput',true),2);
    temp_stim2move_error= std(cellfun(@mean,temp_stim2move,'UniformOutput',true),0,2)./sqrt(size(temp_stim2move,2));
    temp_vel=feval(@(c) cat(2,c{:}), wf_stim_kernels_concat.wheel_velocity(ismember(wf_stim_kernels_concat.group,groups{curr_group})));
    temp_vel2= feval(@(c)  arrayfun(@(id)  cat(2,c{id,:}) ,1:size(temp_vel,1),'uni',false),cellfun(@(x)  median(x,1,'omitmissing')',temp_vel,'UniformOutput',false));
    temp_vel_mean=cellfun(@(x)  median(x,2,'omitmissing'),temp_vel2,'UniformOutput',false  );
    temp_vel_error=cellfun(@(x)  std(x,0,2)./sqrt(size(x,2)),temp_vel2,'UniformOutput',false  );

    % tem_image_kernels_mean=feval(@(c) cat(4,c{:}), arrayfun(@(id)  nanmean(cat(4,tem_image{id,:}),4)         ,1:size(temp_vel,1),'UniformOutput',false))

    tem_trace=cellfun(@(x)   ds.make_each_roi(x,t_kernels,roi1),tem_image,'UniformOutput',false);
    tem_trace_mean=arrayfun(@(id)    median(cat(3,tem_trace{id,:}),3,'omitmissing')  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_error=arrayfun(@(id)    std(cat(3,tem_trace{id,:}),0,3)./sqrt(size(tem_trace,2))  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_peak=arrayfun(@(id)   permute(max(cat(3,tem_trace{id,:}),[],2),[3,1,2])  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_peak_error=arrayfun(@(id)    std(max(cat(3,tem_trace{id,:}),[],2),0,3)./sqrt(size(tem_trace,2))  ,1:size(temp_vel,1),'UniformOutput',false);

    % figure('Position',[50 50 200 600]);




    temp_wf_raw=feval(@(a)  cat(1,a{:}), cellfun(@(a) arrayfun(@(id)  nanmean(cat(3,a{id,:}),3) ,1:4,'uni',false)  ,...
        cellfun(@(x)  cat(2,x{:})         ,...
        wf_stim_kernels_concat.wf_raw(ismember(wf_stim_kernels_concat.group,groups{curr_group})),'uni',false),...
        'UniformOutput',false));

    tem_image_raw=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),temp_wf_raw,'UniformOutput',false);
    tem_trace_raw=cellfun(@(x)   ds.make_each_roi(x,t_kernels,roi1),tem_image_raw,'UniformOutput',false);
    tem_image_raw_mean= arrayfun(@(a) nanmean(cat(3,tem_trace_raw{:,a}),3) ,1:4,'uni',false);
    tem_image_raw_error=arrayfun(@(a) nanstd(cat(3,tem_trace_raw{:,a}),0,3)./sqrt(size(tem_trace_raw,1)) ,1:4,'uni',false);
    tem_image_nanmean=feval(@(a)   cat(4,a{:}), arrayfun(@(a) nanmean(cat(4,tem_image_raw{:,a}),4) ,1:4,'uni',false));



    % tem_image_raw_mean= feval(@(b)   cat(4,b{:}), arrayfun(@(a) nanmean(cat(4,tem_image_raw{:,a}),4) ,1:4,'uni',false));
    % ap.imscroll(tem_image_raw_mean,t_task)
    % axis image off
    % clim( 0.02*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['KWB']));

    % figure('Position',[50 50 200 600]);
    imageLayout= tiledlayout(mainLayout,5,1,'Padding','tight','TileSpacing','none','TileIndexing','rowmajor')
    imageLayout.Layout.Tile = 2*curr_group-1;  % 明确放在主 layout 的第 1 个 tile

    for curr_roi=[used_roi(curr_group) 1 3  14]
        nexttile(imageLayout)
        for curr_state=2:size(temp_vel,1)
            hold on
            ap.errorfill(t_task,tem_image_raw_mean{curr_state}(curr_roi,:),...
                tem_image_raw_error{curr_state}(curr_roi,:),Colors{curr_state},0.5);
            xlim([-0.1 1])
            % if curr_roi

            if curr_roi==7 ||curr_roi==10
                temp_max=feval(@(a) max(a(curr_roi,:,:),[],'all'),cat(5,tem_image_raw_mean{:}))*1.2;
                temp_min=feval(@(a) min(a(curr_roi,:,:),[],'all'),cat(5,tem_image_raw_mean{:}))*1.5;
            else
                temp_max=0.015;
                temp_min=-0.005;
            end

            ylim([temp_min temp_max ])
            % ylim([-0.0002 0.0003])
            xline(0,'.k')
            xline(temp_stim2move_mean(curr_state),'Color',Colors{curr_state})
            axis off
        end


    end

    nexttile(imageLayout)
    hold on
    cellfun(@(x,y,z) ap.errorfill(surround_time_points,x,y,z,0.5),...
        temp_vel_mean(2:4),temp_vel_error(2:4),Colors(2:4)','UniformOutput',false)
    temp_max2 = feval(@(a)  max(a,[],'all') ,  cat(1,temp_vel_mean{:}));
    temp_min2 = feval(@(a)  min(a,[],'all') ,  cat(1,temp_vel_mean{:}))*1.2;
    ylim([temp_min2 temp_max2])
    % ylim([-1500 2000])
    xlim([-0.1 1])
    xline(temp_stim2move_mean(curr_state),'.r')
    xline(0,'.k')
    axis off


    imageLayout=tiledlayout(mainLayout,5,1,'Padding','tight','TileSpacing','none','TileIndexing','rowmajor')

    imageLayout.Layout.Tile = 2*curr_group;  % 明确放在主 layout 的第 1 个 tile

    for curr_roi=[used_roi(curr_group) 1 3  14]
        nexttile(imageLayout)
        for curr_state=2:size(temp_vel,1)
            hold on
            ap.errorfill(t_kernels,tem_trace_mean{curr_state}(curr_roi,:),...
                tem_trace_error{curr_state}(curr_roi,:),Colors{curr_state},0.1);
            xlim([-0.1 1])
            % if curr_roi

            if curr_roi==7 ||curr_roi==10
                temp_max=feval(@(a) max(a(curr_roi,:,:),[],'all'),cat(5,tem_trace_mean{:}))*1.2;
                temp_min=feval(@(a) min(a(curr_roi,:,:),[],'all'),cat(5,tem_trace_mean{:}))*1.5;
            else
                temp_max=0.0004;
                temp_min=-0.0001;
            end

            ylim([temp_min temp_max ])
            % ylim([-0.0002 0.0003])
            xline(0,'.k')
            xline(temp_stim2move_mean(curr_state),'Color',Colors{curr_state})
            axis off
        end
        if curr_group==2 &curr_roi==used_roi(curr_group)
            h = findobj(gca,'Type','line');
            legend(h([1 2 3]), {'>0.3s','0.2-0.3s','0-0.2s'},'Box','off')
        end
    end

    nexttile(imageLayout)
    hold on
    cellfun(@(x,y,z) ap.errorfill(surround_time_points,x,y,z,0.5),...
        temp_vel_mean(2:4),temp_vel_error(2:4),Colors(2:4)','UniformOutput',false)
    temp_max2 = feval(@(a)  max(a,[],'all') ,  cat(1,temp_vel_mean{:}));
    temp_min2 = feval(@(a)  min(a,[],'all') ,  cat(1,temp_vel_mean{:}))*1.2;
    ylim([temp_min2 temp_max2])
    % ylim([-1500 2000])
    xlim([-0.1 1])
    xline(temp_stim2move_mean(curr_state),'.r')
    xline(0,'.k')
    axis off




end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_A2.eps'), ...
    'ContentType','vector');

%% EDF 3
clear all
clc
Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';

% Path='D:\Data process\slide\papers';
U_master = plab.wf.load_master_U;
% load('C:\Users\dsong\Documents\MATLAB\Da_Song\DS_scripts_ptereslab\General_information\roi.mat')
load(fullfile(Path,'data\General_information\roi.mat'))
surround_samplerate = 35;
surround_window_task = [-0.2,1];
task_boundary1=0;
task_boundary2=0.2;

t_kernels=1/surround_samplerate*[-10:30];
kernels_period=find(t_kernels>task_boundary1&t_kernels<task_boundary2);



main_preload_vars = who;
load(fullfile(Path,'data','wf_task_kernels'));
tem_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[7 8],:)),  wf_task_kernels_across_day,'UniformOutput',false);
image_acorss_time=cellfun(@(x)    nanmean(x,[4 5]),tem_image,'UniformOutput',false  );

% images_max=cellfun(@(x)   nanmean(max(x(:,:,kernels_period,:,:),[],3),[4 5]   ),tem_image,'UniformOutput',false  );
images_max=cellfun(@(x)  max(x(:,:,kernels_period),[],3),image_acorss_time,'UniformOutput',false  );


scale_image=0.0003;
Color={'B','R'};

figure('Position',[50 50 900 300])
tiledlayout (3,length(find(t_kernels>-0.05& t_kernels<0.25)),'TileSpacing','none')

for curr_group=1:2

    for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
        ax=nexttile
        imagesc(image_acorss_time{curr_group}(:,:,curr_frame))
        axis image off;
        clim(scale_image .* [0, 1]);
        colormap(ax, ap.colormap(['W' Color{curr_group}] ));
        ap.wf_draw('ccf', [0.5 0.5 0.5]);
        if curr_group==1
            title( num2str(t_kernels(curr_frame),'%.2f'),'FontWeight','normal')
        end

    end
end

for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
    ax=nexttile

    image(ds.colormap_overlay(image_acorss_time{1}(:,:,curr_frame),...
        image_acorss_time{2}(:,:,curr_frame),...
        'B','R',[],scale_image))

    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end

exportgraphics(gcf, fullfile(Path,'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig EDF3_a.eps'), ...
    'ContentType','vector');

figure('Position', [50 50 200 200] );
image(ds.colormap_overlay(images_max{1},...
    images_max{2},...
    'B','R',[],scale_image))
axis image off;
ap.wf_draw('ccf', [0.5 0.5 0.5]);
exportgraphics(gcf, fullfile(Path,'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig EDF3_b.eps'), ...
    'ContentType','vector');

load(fullfile(Path,'data','wf_passive_kernels'));

tem_image_pass=cellfun(@(x,id) plab.wf.svd2px(U_master(:,:,1:size(x{id},1)),x{id}(:,:,:,[10 11],:)), ...
    wf_passive_kernels_across_day,{1;2},'UniformOutput',false);

image_pass_acorss_time=cellfun(@(x)    nanmean(x,[5 6]),tem_image_pass,'UniformOutput',false  );
% iamges_passive_max=cellfun(@(x)    nanmean(max(x(:,:,kernels_period,:,:,:),[],[3,4]),[5 6]),tem_image_pass,'UniformOutput',false  );
iamges_passive_max=cellfun(@(x)    max(x(:,:,kernels_period,:),[],[3,4]),image_pass_acorss_time,'UniformOutput',false  );

figure('Position', [50 50 900 700] )
mainfig=tiledlayout(7,1,'TileSpacing','none')
for curr_group=1:2
    subfig=tiledlayout(mainfig,3,sum(t_kernels>-0.05& t_kernels<0.25),'TileSpacing','none');
    subfig.Layout.Tile=3*curr_group-2;
    subfig.Layout.TileSpan = [3 1];      % 跨 2 行 × 1 列

    for curr_stim=1:3
        for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
            ax=nexttile(subfig);
            imagesc(image_pass_acorss_time{curr_group}(:,:,curr_frame,curr_stim));
            axis image off;
            clim(scale_image .* [0, 1]);
            colormap(ax, ap.colormap(['W' Color{curr_group}] ));
            clim(scale_image .* [0, 1]);
            ap.wf_draw('ccf', [0.5 0.5 0.5]);

        end
    end
end
subfig=tiledlayout(mainfig,1,sum(t_kernels>-0.05& t_kernels<0.25),'TileSpacing','none');
subfig.Layout.Tile=7;

for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
    ax=nexttile(subfig);

    image(ds.colormap_overlay(image_pass_acorss_time{1}(:,:,curr_frame,3),...
        image_pass_acorss_time{2}(:,:,curr_frame,2),...
        'B','R',[],scale_image))

    axis image off;

    ap.wf_draw('ccf', [0.5 0.5 0.5]);

end

exportgraphics(gcf, fullfile(Path,'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig EDF3_c.eps'), ...
    'ContentType','vector');
figure('Position', [50 50 200 200] );
image(ds.colormap_overlay(iamges_passive_max{1},...
    iamges_passive_max{2},...
    'B','R',[],scale_image))
axis image off;
ap.wf_draw('ccf', [0.5 0.5 0.5]);

exportgraphics(gcf, fullfile(Path,'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig EDF3_d.eps'), ...
    'ContentType','vector');

figure('Position', [50 50 80 80] );
ds.bi_colorbar('B','R',scale_image)
exportgraphics(gcf, fullfile(Path,'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig EDF3_e.eps'), ...
    'ContentType','vector');

%% EDF C correlation between task vs passive kernels
clear all
Path='D:\Data process\project_cross_model\wf_data\data_package';

U_master = plab.wf.load_master_U;
load('C:\Users\dsong\Documents\MATLAB\Da_Song\DS_scripts_ptereslab\General_information\roi.mat');
temp_data_all=cell(2,1);
surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
t_kernels=[-10:30]/surround_samplerate;
period=find(t_kernels>0&t_kernels<0.2);
surround_window = [-0.5,1];
mousecam_framerate = 30;
face_time = surround_window(1):1/mousecam_framerate:surround_window(2);

for curr_group=1:2
    switch curr_group
        case 1
            animals = {'DS007','DS010','AP019','AP021','DS011','AP022'};
        case 2
            animals = {'DS000','DS004','DS014','DS015','DS016'};
    end

    temp_data_all{curr_group}=table;
    for curr_animal=1:length(animals)
        preload_vars=who;
        animal=animals{curr_animal};
        data_all=load(fullfile(Path,[animal '_all_data.mat']));

        switch curr_group
            case 1
                select_id_3=(strcmp([data_all.task_name],'stim_wheel_right_stage1')|strcmp([data_all.task_name],'stim_wheel_right_stage2'))&...
                    ~cellfun(@isempty ,data_all.wf_lcr_passive);
            case 2
                select_id_3=(strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume')|strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume'))&...
                    ~cellfun(@isempty ,data_all.wf_hml_passive_audio);


        end
        temp_p_val=arrayfun(@(id)  data_all.behavior_task{id}.rxn_l_p(1)<0.05, find(select_id_3),'UniformOutput',true);


        %% wf_task

        all_wf_task=[data_all.wf_task(select_id_3)];
        % all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.iti_move_kernels{1}),all_wf_task,'UniformOutput',false);
        all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.stim_kernels{1}),all_wf_task,'UniformOutput',false);
        % all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.all_iti_move_kernels{1}),all_wf_task,'UniformOutput',false);

        temp_data_all{curr_group}.task_image{curr_animal}=cat(4,all_task_image{find(temp_p_val==1,2,'last')});
        %% wf_passive
        switch curr_group
            case 1
                all_wf_passive_1=[data_all.wf_lcr_passive(select_id_3)];
                used_passive=3;
            case 2
                all_wf_passive_1=[data_all.wf_hml_passive_audio(select_id_3)];
                used_passive=2;
        end
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.kernels_decoding{used_passive},1)),x.kernels_decoding{used_passive}),all_wf_passive_1,'UniformOutput',false);

        temp_data_all{curr_group}.passive_image{curr_animal}=cat(4,all_passive_image_pre{find(temp_p_val==1,2,'last')})

        clearvars('-except',preload_vars{:});
    end

end

pairs={[1 1 ],[2 2],[1 2],[2 1]}
pairs_names={'Visual task vs Visual passive','Auditory task vs Auditory passive',...
    'Visual task vs Auditory passive','Auditory task vs Visual passive'};
temp_val=cell(4,1);
figure('Position',[50 50 600 900])
for curr_pair=1:length(pairs)


    % task_images_max=feval(@(a)  cat(3,a{:}), ...
    %     cellfun(@(x)  permute(max(x(:,:,period,:),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(1)}.task_image,'UniformOutput',false));
    %
    task_images_max=feval(@(a)  cat(3,a{:}), ...
        cellfun(@(x)  permute(max(nanmean(x(:,:,period,:),4),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(1)}.task_image,'UniformOutput',false));


    task_length=size(task_images_max,3);

    % passive_images_max=feval(@(a)  cat(3,a{:}), ...
    %     cellfun(@(x)  permute(max(x(:,:,period,:),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(2)}.passive_image,'UniformOutput',false));
    %

    passive_images_max=feval(@(a)  cat(3,a{:}), ...
        cellfun(@(x)  permute(max(nanmean(x(:,:,period,:),4),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(2)}.passive_image,'UniformOutput',false));

    passive_length=size(passive_images_max,3);


    X=cat(3,task_images_max,passive_images_max);


    % ap.imscroll(X)
    % axis image off;
    % ap.wf_draw('ccf', [0.5 0.5 0.5]);
    % clim(0.0003 .* [ 0, 1]);
    % colormap( ap.colormap('WB' ));


    X = double(X);   % n*m*10

    K = size(X, 3);  % 这里应该是 10
    C = zeros(K, K);

    for i = 1:K
        vi = X(:,:,i);
        vi = vi(:);

        for j = 1:K
            vj = X(:,:,j);
            vj = vj(:);

            % 去掉常数向量导致的 NaN
            if std(vi) == 0 || std(vj) == 0
                C(i,j) = NaN;
            else
                C(i,j) = corr(vi, vj);
            end
        end
    end

    temp_val{curr_pair}=C(1:task_length,task_length+1:task_length+passive_length);


    K = size(C,1);

    % 创建mask：保留对角线 + 上三角
    mask = triu(true(K));

    % 其余设为 NaN（不显示）
    C_plot = nan(K);
    C_plot(mask) = C(mask);


    nexttile
    imagesc(C_plot);
    colormap( ap.colormap('WB' ));
    set(gca, 'XAxisLocation', 'top');
    set(gca, 'YAxisLocation', 'right');
    axis image;
    box off

    labels={'task','passive'}
    set(gca, 'XTick', [task_length/2 task_length+passive_length/2], 'YTick', [task_length/2 task_length+passive_length/2]);
    set(gca, 'XTickLabel', labels, 'YTickLabel', labels);
    % set(gca, 'TickLabelInterpreter', 'none');

    line([task_length+0.5 (task_length+passive_length)+0.5 ], [task_length+0.5 task_length+0.5 ],'Color',[1 0 0]);
    line([ task_length+0.5 task_length+0.5],[0.5 task_length+0.5 ],'Color',[1 0 0]);
    title(pairs_names{curr_pair})
end
colorbar;
caxis([0 1]);   % correlation 范围

temp_v=cellfun(@(x)  x(:),temp_val,'UniformOutput',false)
% figure('Position',[50 50 300 400])
nexttile
ds.make_bar_plot(cellfun(@(x)  x(:),temp_val,'UniformOutput',false),'ShowDots',0)
set(gca, 'XTick',[1:4],'XTickLabel', pairs_names);
ylim([0 1])
ylabel('Correlation')

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_C1.eps'), ...
    'ContentType','vector');

ds.shuffle_test(temp_v{3},temp_v{4})



% correlation trace
image_corr=cell(2,1);
for curr_group=1:2
    task_images=cat(4,temp_data_all{curr_group}.task_image{:});
    passive_images=cat(4,temp_data_all{curr_group}.passive_image{:});
    % image_corr=cell(size(task_images,4),1);
    for curr_image=1:size(task_images,4)
        used_task=reshape(task_images(:,:,:,curr_image),[],size(task_images,3))';
        used_passive=reshape(passive_images(:,:,:,curr_image),[],size(passive_images,3))';
        Cmat=zeros(size(used_task,2),1);
        for curr_pixel=1:size(used_task,2)
            Cmat(curr_pixel) = corr(used_task(:,curr_pixel),used_passive(:,curr_pixel));   % 注意转置！
        end
        image_corr{curr_group}{curr_image} = reshape(Cmat, size(task_images,1), size(task_images,2));
    end
end

figure('Position',[50 50 400 200]);
for curr_group=1:2
    nexttile
    plot_image=nanmean(cat(3,image_corr{curr_group}{:}),3);
    % plot_image(plot_image<0.3)=0;
    imagesc(plot_image)
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim([ 0, 1]);
    colormap( ap.colormap('WB' ));
end

colorbar
exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_C2.eps'), ...
    'ContentType','vector');

%% EDF D audio balance to visual


animals = {'DS029','DS030','DS031'};

temp_data_all=table;
for curr_animal=1:length(animals)
    preload_vars=who;

    animal=animals{curr_animal};
    data_all=load(fullfile(Path,[animal '_all_data.mat']));


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
all_audio_passive_image_max_pre=permute(max(nanmean(all_passive_image_pre(:,:,period,:,:),5),[],3),[1,2,4,3]);
temp_wf_passive_plot_tace= ds.make_each_roi(all_passive_image_pre, length(t_kernels),roi1);


temp_hml_passive_post= feval(@(a) cat(4,a{:}) ,cellfun(@(x) nanmean(cat(4,x{:}),4),temp_data_all.all_hml_passive_post,'UniformOutput',false  ));
all_passive_image_post= plab.wf.svd2px(U_master(:,:,1:size(cat(3,temp_hml_passive_post),1)),temp_hml_passive_post);
all_audio_passive_image_max_post=permute(max(nanmean(all_passive_image_post(:,:,period,:,:),5),[],3),[1,2,4,3]);
temp_wf_passive_plot_tace= ds.make_each_roi(all_passive_image_pre, length(t_kernels),roi1);


temp_lcr_passive= feval(@(a) cat(4,a{:}) ,cellfun(@(x) nanmean(cat(4,x{:}),4),temp_data_all.all_lcr_passive,'UniformOutput',false  ));
all_visual_passive_image= plab.wf.svd2px(U_master(:,:,1:size(cat(3,temp_lcr_passive),1)),temp_lcr_passive);
all_visual_passive_image_max=permute(max(nanmean(all_visual_passive_image(:,:,period,:),5),[],3),[1,2,4,3]);
temp_visual_passive_plot_tace= ds.make_each_roi(all_visual_passive_image, length(t_kernels),roi1);

temp_visual=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_visual_task{:})));
all_visual_task_image= plab.wf.svd2px(U_master(:,:,1:size(temp_visual,1)),temp_visual);
all_visual_task_image_max=max(nanmean(all_visual_task_image(:,:,period,:),4),[],3);

temp_audio_pre=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_audio_task_pre{:})));
all_audio_task_image_pre= plab.wf.svd2px(U_master(:,:,1:size(temp_audio_pre,1)),temp_audio_pre);
all_audio_task_image_max_pre=max(nanmean(all_audio_task_image_pre(:,:,period,:),4),[],3);


temp_audio_post=feval(@(b)   cat(3,b{:})  ,feval(@(a)  cat(1,a{:}) , cat(1,temp_data_all.all_audio_task_post{:})));
all_audio_task_image_post= plab.wf.svd2px(U_master(:,:,1:size(temp_audio_post,1)),temp_audio_post);
all_audio_task_image_max_post=max(nanmean(all_audio_task_image_post(:,:,period,:),4),[],3);


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
%%  EDF E  different visual task types
clear all
clc
U_master = plab.wf.load_master_U;
load(fullfile('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data\General_information\roi.mat'))


load(fullfile('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data\revision',...
    'different_task_behavior.mat'))

surround_samplerate = 35;
t_kernels=1/surround_samplerate*[-10:30];
%

colors={[0.4 0.4 1],[0.8 0.8 1],[0.5 0.5 1],[0.44 0.4 1],[0.44 0.4 1]};
use_group={[1 2],3}
for curr_groups=1:2
    figure('Position',[50 50 600 150]);
    tiledlayout(1,4)
    nexttile
    hold on
    for curr_group=use_group{curr_groups}
        cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.reaction_time{curr_group},'UniformOutput',false )

        cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.reaction_time{curr_group},...
            behavior_data.p_val{curr_group},'UniformOutput',false )
    end
    % ap.errorfill(1:length(RT_mean),RT_mean,RT_error,[0 0 0])
    xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
    yticks([0.1 0.5 1 2 4 6])
    set(gca, 'YScale', 'log')
    ylim([0.1 10])

    ylabel('Reaction time (s)')
    xlabel('Days')
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    nexttile
    hold on
    for curr_group=use_group{curr_groups}

        cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.performance{curr_group},'UniformOutput',false )

        cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.performance{curr_group},...
            behavior_data.p_val{curr_group},'UniformOutput',false )
    end

    xlim([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])
    ylim([-0.2 1])
    ylabel('Performance')
    xlabel('Days')
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])
    yticks([0 1])

    % sgtitle(title_name)


    temp_wf_task=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(3,a{end-1:end}),3)...
        ,x,'UniformOutput',false) ,behavior_data.task_kernels(use_group{curr_groups}),'UniformOutput',false)

    tem_image_task=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
        ,temp_wf_task,'UniformOutput',false);

    image_plot_task=feval(@(w)    nanmean(max(w(:,:,t_kernels>0&t_kernels<0.2 ,:),[],3),4) ,...
        feval(@(c) cat(4,c{:}), cat(2,tem_image_task{:})));

    nexttile
    imagesc(image_plot_task)
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));



    temp_wf_passive=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(4,a{end-1:end}),4)...
        ,x,'UniformOutput',false) ,behavior_data.passive_kernels(use_group{curr_groups}),'UniformOutput',false)

    tem_image_passive=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1),:),x),a,'UniformOutput',false)...
        ,temp_wf_passive,'UniformOutput',false);

    image_plot_passive=feval(@(w)    permute(nanmean(max(w(:,:,t_kernels>0&t_kernels<0.2 ,:,:),[],3),5),[1,2,4,3]) ,...
        feval(@(c) cat(5,c{:}), cat(2,tem_image_passive{:})));

    nexttile
    imagesc(image_plot_passive(:,:,3))
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));

    exportgraphics(gcf, fullfile(plab.locations.server_path,...
        ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_E' num2str(curr_groups) '.eps']), ...
        'ContentType','vector');
end

use_group=[ 4 5];
figure('Position',[50 50 1000 150]);
mainfig=tiledlayout(1,4)

plot_fig=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'none', 'Padding', 'none');
plot_fig.Layout.Tile = 1;  % 明确放在主 layout 的第 1 个 tile

for curr_group=use_group
    nexttile(plot_fig)
    hold on
    cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.reaction_time{curr_group},'UniformOutput',false )

    cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.reaction_time{curr_group},...
        behavior_data.p_val{curr_group},'UniformOutput',false )
    % ap.errorfill(1:length(RT_mean),RT_mean,RT_error,[0 0 0])
    set(gca, 'YScale', 'log')

    ylim([0.1 10])
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    switch curr_group
        case 4
            xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            yticks([0.1 0.5 1 5 10])
            ylabel('Reaction time (s)')
            xlabel('Days')
        case 5
            xlim([0 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            set(gca,'YColor','none')
    end
end

plot_fig2=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'none', 'Padding', 'none');
plot_fig2.Layout.Tile = 2;  % 明确放在主 layout 的第 1 个 tile
for curr_group=use_group
    nexttile (plot_fig2)
    hold on
    cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.performance{curr_group},'UniformOutput',false )
    cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.performance{curr_group},...
        behavior_data.p_val{curr_group},'UniformOutput',false )

    ylim([-0.2 1])
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    switch curr_group
        case 4
            xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            yticks([0.1 0.5 1 2 4 6])
            % set(gca, 'YScale', 'log')
            ylabel('Performance')
            xlabel('Days')
        case 5
            xlim([0 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            set(gca,'YColor','none')
    end
end
% sgtitle(title_name)


temp_wf_task=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(3,a{end-1:end}),3)...
    ,x,'UniformOutput',false) ,behavior_data.task_kernels(use_group),'UniformOutput',false)
tem_image_task=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
    ,temp_wf_task,'UniformOutput',false);
image_plot_task=cellfun(@(x) feval(@(v)  nanmean(v,3) , feval(@(c) cat(3,c{:}) ,cellfun(@(a)  max(a(:,:,t_kernels>0&t_kernels<0.2 ,:),[],3),x,'UniformOutput',false)))...
    ,tem_image_task,'UniformOutput',false  );

plot_fig3=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'tight', 'Padding', 'tight');
plot_fig3.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile

for curr_image=1:2
    nexttile (plot_fig3)
    imagesc(image_plot_task{curr_image})
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));
end


temp_wf_passive=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(4,a{end-1:end}),4)...
    ,x,'UniformOutput',false) ,behavior_data.passive_kernels(use_group),'UniformOutput',false)
tem_image_passive=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1),:),x),a,'UniformOutput',false)...
    ,temp_wf_passive,'UniformOutput',false);
image_plot_passive=cellfun(@(x) feval(@(v)  nanmean(v,4) , feval(@(c) cat(4,c{:}) ,cellfun(@(a)  permute (max(a(:,:,t_kernels>0&t_kernels<0.2 ,:,:),[],3),[1,2,4,3]),x,'UniformOutput',false)))...
    ,tem_image_passive,'UniformOutput',false  );
plot_fig4=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'tight', 'Padding', 'tight');
plot_fig4.Layout.Tile = 4;  % 明确放在主 layout 的第 1 个 tile

for curr_image=1:2
    nexttile (plot_fig4)
    imagesc(image_plot_passive{curr_image}(:,:,3))
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));
end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_E3.eps'), ...
    'ContentType','vector');

%% EDF F   viarable task  contrast volume
clear all
U_master=plab.wf.load_master_U;
load(fullfile('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data\General_information\roi.mat'))

tempdata{1}= load(fullfile(plab.locations.server_path,'Lab\Papers\Song_2025\data\revision\visual_task_variable.mat'));
tempdata{2}=load(fullfile(plab.locations.server_path,'Lab\Papers\Song_2025\data\revision\audio_task_variable.mat'));

surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
t_kernels=[-10:30]/surround_samplerate;
period=find(t_kernels>0&t_kernels<0.2);


temp_perform=cellfun(@(s) feval(@(a)  cat(2,a{:}) ,cellfun(@(x) ...
    nanmean(x,2), s.kernels_data.performance,'UniformOutput',false)),tempdata,'uni',false     );

temp_perform_mean=cellfun(@(x) mean(x,2),temp_perform,'UniformOutput',false);
temp_perform_error=cellfun(@(x) std(x,0,2)./sqrt(size(x,2)),temp_perform,'UniformOutput',false);



temp_rxt=cellfun(@(s) feval(@(a)  cat(2,a{:}) ,cellfun(@(x) ...
    nanmean(x,2), s.kernels_data.react_time,'UniformOutput',false)),tempdata,'uni',false     );

temp_rxt_mean=cellfun(@(x) mean(x,2),temp_rxt,'UniformOutput',false);
temp_rxt_error=cellfun(@(x) std(x,0,2)./sqrt(size(x,2)),temp_rxt,'UniformOutput',false);



figure('Position',[50 50 400 200])
nexttile
hold on
ap.errorfill(1:6,temp_rxt_mean{1},temp_rxt_error{1},[0 0 1],0.1);
ap.errorfill(1:6,temp_rxt_mean{2},temp_rxt_error{2},[1 0 0],0.1);
xticks([1 6]);
xticklabels({'min','max'}); % 先清掉默认标签

xlabel('Contrast/Volume')
xlim([1 6])
ylim([0 6])
ylabel('Reaction time (s)');
set(gca,'Color','none', 'YScale', 'log')

nexttile
hold on
ap.errorfill(1:6,temp_perform_mean{1},temp_perform_error{1},[0 0 1],0.1);
ap.errorfill(1:6,temp_perform_mean{2},temp_perform_error{2},[1 0 0],0.1);
xticks([1 6]);
xticklabels({'min','max'}); % 先清掉默认标签

xlabel('Contrast/Volume')
xlim([1 6])
ylim([0 1])
ylabel('Performance');
set(gca,'Color','none')

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_F1.eps'), ...
    'ContentType','vector');
%

color_name={'B','R'};
figure('Position',[50 50 600 200])
tiledlayout(2,6,'TileSpacing','none')
for curr_group=1:2
    kernels_data=tempdata{curr_group}.kernels_data;
    temp_images=feval(@(a) cat(5,a{:}),...
        cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),kernels_data.kernels,'UniformOutput',false));
    temp_image_max=permute(nanmean(max(temp_images(:,:,period,:,:),[],3),5),[1,2,4,3,5]);
    for curr_state=1:size(temp_image_max,3)
        a1=nexttile
        imagesc(temp_image_max(:,:,curr_state))
        clim( 0.0003*[0,1]);
        ap.wf_draw('ccf',[0.5 0.5 0.5]);
        colormap( a1,ap.colormap(['W' color_name{curr_group} ]));
        axis image off
    end
end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_F2.eps'), ...
    'ContentType','vector');
%
figure('Position',[50 50 400 300])
nexttile
sensory_roi=[7 9]
cmap{1}=colormap([linspace(0.7,0,6)', linspace(0.7,0,6)', linspace(1,1,6)']) % 浅蓝→深蓝
cmap{2}=colormap([linspace(1,1,6)', linspace(0.7,0,6)', linspace(0.7,0,6)'])   % 浅红→深红
R_task=[];
P_task=[];
for curr_group=1:2

    kernels_data=tempdata{curr_group}.kernels_data;
    temp_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),kernels_data.kernels,'UniformOutput',false);
    temp_plot_tace=feval(@(a) cat(4,a{:}), cellfun(@(x)  ds.make_each_roi( nanmean(x,5), length(t_kernels),roi1)  , temp_image,'UniformOutput',false  ));
    temp_plot_tace_max=permute(max(temp_plot_tace(:,period,:,:),[],2),[3,1,4,2]);

    temp_x=nanmean(temp_plot_tace_max(:,sensory_roi(curr_group),:),3);
    temp_y= nanmean(temp_plot_tace_max(:,1,:),3);
    hold on
    scatter(temp_x, temp_y, 50, cmap{curr_group},'filled')


    p_passive = polyfit(temp_x, temp_y, 1);
    x_fit_passive =[0 0.0002] ;
    y_fit_passive = polyval(p_passive, x_fit_passive);
    % plot(x_fit_passive, y_fit_passive, '-', 'LineWidth', 2,'Color',cmap{curr_group}(1,:));
    [R_task(curr_group),P_task(curr_group)] = corr(temp_x,temp_y);
    box off




end

ylabel('mPFC \DeltaF/F_0')
xlabel( ' Sensory areas \DeltaF/F_0')
axis square
set(gca,'Color','none')

nexttile
% figure('Position',[50 50 200 300])
sensory_roi=[7 9]
cmap{1}=colormap([linspace(0.7,0,6)', linspace(0.7,0,6)', linspace(1,1,6)']) % 浅蓝→深蓝
cmap{2}=colormap([linspace(1,1,6)', linspace(0.7,0,6)', linspace(0.7,0,6)'])   % 浅红→深红
R_task=[];
P_task=[];
for curr_group=1:2

    kernels_data=tempdata{curr_group}.kernels_data;
    temp_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),kernels_data.kernels,'UniformOutput',false);
    temp_plot_tace=feval(@(a) cat(4,a{:}), cellfun(@(x)  ds.make_each_roi( nanmean(x,5), length(t_kernels),roi1)  , temp_image,'UniformOutput',false  ));
    temp_plot_tace_max=permute(max(temp_plot_tace(:,period,:,:),[],2),[3,1,4,2]);

    temp_x=nanmean(temp_plot_tace_max(:,sensory_roi(curr_group),:),3);
    temp_y= nanmean(temp_plot_tace_max(:,3,:),3);
    hold on
    scatter(temp_x, temp_y, 50, cmap{curr_group},'filled')


    p_passive = polyfit(temp_x, temp_y, 1);
    x_fit_passive =[0 0.0002] ;
    y_fit_passive = polyval(p_passive, x_fit_passive);
    % plot(x_fit_passive, y_fit_passive, '-', 'LineWidth', 2,'Color',cmap{curr_group}(1,:));
    [R_task(curr_group),P_task(curr_group)] = corr(temp_x,temp_y);
    box off


end
ylabel('aPFC \DeltaF/F_0')
xlabel( ' Sensory areas \DeltaF/F_0')
axis square
set(gca,'Color','none')

% ===== 自定义 legend（右上角外侧）=====
ax_leg = axes('Position',[0.2 0.8 0.3 0.18]); % ← 关键：>0.8 基本就在外面了
hold(ax_leg,'on')
axis(ax_leg,'off')

% data1（蓝色）
for i = 1:6
    scatter(ax_leg, i, 2, 40, cmap{1}(i,:), 'filled')
end
text(ax_leg, 7, 2, 'Visual contrast', 'VerticalAlignment','middle')

% data2（红色）
for i = 1:6
    scatter(ax_leg, i, 1, 40, cmap{2}(i,:), 'filled')
end
text(ax_leg, 7, 1, 'Auditory volume', 'VerticalAlignment','middle')

xlim(ax_leg,[0 8])
ylim(ax_leg,[0.5 2.5])

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_F3.eps'), ...
    'ContentType','vector');

%% EDF G  visual passive size 20 60

clear all

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025\data';
load(fullfile(save_path,'revision','visual_size_passive_compare.mat'));
surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
t_kernels=[-10:30]/surround_samplerate;
period=find(t_kernels>0&t_kernels<0.2);
load('C:\Users\dsong\Documents\MATLAB\Da_Song\DS_scripts_ptereslab\General_information\roi.mat')

U_master = plab.wf.load_master_U;


tem_passive_s_image=cell(2,1);
tem_passive_l_image=cell(2,1);
trace_s_mean=cell(2,1);
trace_l_mean=cell(2,1);
for curr_animal=1:5

    tem_passive_s=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),passive_data.lcr_passive(curr_animal),'UniformOutput',false);
    tem_passive_l=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),passive_data.lcr_passive_size60(curr_animal),'UniformOutput',false);
    tem_passive_s_image{curr_animal}=permute(nanmean(max(tem_passive_s{1}(:,:,period,3,:),[],3),5),[1,2,5,3,4]);
    tem_passive_l_image{curr_animal}=permute(nanmean(max(tem_passive_l{1}(:,:,period,3,:),[],3),5),[1,2,5,3,4]);

    trace_s=ds.make_each_roi(permute(tem_passive_s{1}(:,:,:,3,:),[1,2,3,5,4]),t_kernels,roi1);
    trace_l=ds.make_each_roi(permute(tem_passive_l{1}(:,:,:,3,:),[1,2,3,5,4]),t_kernels,roi1);
    trace_s_mean{curr_animal}=nanmean(trace_s,3);
    trace_s_error=std(trace_s,0,3,'omitmissing')./sqrt(size(trace_s,3));
    trace_l_mean{curr_animal}=nanmean(trace_l,3);
    trace_l_error=std(trace_l,0,3,'omitmissing')./sqrt(size(trace_l,3));

    % 
    % figure('Position',[50 50 800 200]);
    % tiledlayout(1,4)
    % 
    % nexttile
    % imagesc(tem_passive_s_image{curr_animal})
    % axis image off
    % clim( 0.0003*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['PWG']));
    % nexttile
    % imagesc(tem_passive_l_image{curr_animal})
    % axis image off
    % clim( 0.0003*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['PWG']));
    % 
    % nexttile
    % hold on
    % ap.errorfill(t_kernels,trace_s_mean{curr_animal}(7,:),trace_s_error(1,:))
    % ap.errorfill(t_kernels,trace_l_mean{curr_animal}(7,:),trace_l_error(1,:))
    % title('V1')
    % xlim([-0.1 0.5])
    % 
    % nexttile
    % hold on
    % ap.errorfill(t_kernels,trace_s_mean{curr_animal}(1,:),trace_s_error(1,:))
    % ap.errorfill(t_kernels,trace_l_mean{curr_animal}(1,:),trace_l_error(1,:))
    % title('mPFC')
    % xlim([-0.1 0.5])


end

trace_s_all_mean=nanmean(cat(3,trace_s_mean{:}),3);
trace_s_all_error=std(cat(3,trace_s_mean{:}),0,3)./sqrt(length(trace_s_mean));
trace_l_all_mean=nanmean(cat(3,trace_l_mean{:}),3);
trace_l_all_error=std(cat(3,trace_l_mean{:}),0,3)./sqrt(length(trace_l_mean));

tem_passive_l_image_all=nanmean(cat(3,tem_passive_l_image{:}),3);
tem_passive_s_image_all=nanmean(cat(3,tem_passive_s_image{:}),3);



max_data_s=feval(@(a) cat(2,a{:}), cellfun(@(x) max(x([1 7],period),[],2),   trace_s_mean,'UniformOutput',false))
max_data_l=feval(@(a) cat(2,a{:}), cellfun(@(x) max(x([1 7],period),[],2),   trace_l_mean,'UniformOutput',false))

% figure
% hold on
% scatter(max_data_s(2,:), max_data_s(1,:), 50,'filled')
% scatter(max_data_l(2,:), max_data_l(1,:), 50,'filled')
% xlim([0 0.001])
% ylim([0 0.001])

figure('Position',[50 50 800 200]);
tiledlayout(1,4)

nexttile
imagesc(tem_passive_s_image_all)
axis image off
clim( 0.0003*[-1,1]);
ap.wf_draw('ccf',[0.5 0.5 0.5]);
title('size 20')

colormap( ap.colormap(['PWG']));
nexttile
imagesc(tem_passive_l_image_all)
axis image off
clim( 0.0003*[-1,1]);
ap.wf_draw('ccf',[0.5 0.5 0.5]);
colormap( ap.colormap(['PWG']));
title('size 60')

nexttile
hold on
ap.errorfill(t_kernels,trace_s_all_mean(7,:),trace_s_all_error(7,:))
ap.errorfill(t_kernels,trace_l_all_mean(7,:),trace_l_all_error(7,:))
title('V1')
xlim([-0.1 0.5])
set(gca,'Color','none')

nexttile
hold on
ap.errorfill(t_kernels,trace_s_all_mean(1,:),trace_s_all_error(1,:))
ap.errorfill(t_kernels,trace_l_all_mean(1,:),trace_l_all_error(1,:))
title('mPFC')
xlim([-0.1 0.5])
legend({'','Small','','Large'},'Box','off','Location','northeastoutside')
set(gca,'Color','none')


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_G.eps'))

%%     EDF H  block test



load(fullfile(plab.locations.server_path,'Lab\Papers\Song_2025\data\revision\wf_block_test.mat') )

U_master = plab.wf.load_master_U;
load('C:\Users\dsong\Documents\MATLAB\Da_Song\DS_scripts_ptereslab\General_information\roi.mat');

animals =     { 'AP030','AP032','DS030','DS031','DS029'};

temp_image_max=cell(length(animals),1);
temp_plot_tace=cell(length(animals),1);
temp_image_trace=cell(length(animals),1);
for curr_animal=1:length(animals)

    temp_image{1}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task1(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{3}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task2(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{2}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive1(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{4}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive2(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));

    temp_image{5}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task_signle(end),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{6}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive_signle(end),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));

    temp_image_max{curr_animal}= feval(@(x)permute(max(x(:,:,t_kernels>0 & t_kernels<0.2,:),[],3),[1,2,4,3]),cat(4,temp_image{:}));
    temp_plot_tace{curr_animal}= ds.make_each_roi(cat(4,temp_image{:}), length(t_kernels),roi1);

    temp_image_trace{curr_animal}=temp_image;
end


for curr_animal=[1 2 3 4 5]
    ap.imscroll(cat(4,temp_image_trace{curr_animal}{:}),t_kernels)
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim(0.0003 .* [ 0, 1]);
    colormap( ap.colormap('WG' ));
end

figure('Position',[50 50 800 200]);
tiledlayout(1,5)
for curr_image=1:4
    nexttile
    imagesc( feval(@(c) c(:,:,curr_image),  nanmean(cat(4,temp_image_max{:}),4) ));
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim(0.0003 .* [ 0, 1]);
    colormap( ap.colormap('WG' ));


end
plot_trace_mean = feval(@(c) nanmean(c,4 ),cat(4,temp_plot_tace{:}));
plot_trace_error = feval(@(c) std(c,0,4,'omitmissing' )./ sqrt(size(c,4)),cat(4,temp_plot_tace{:}));

plot_max_mean=feval(@(c) permute(nanmean(max(c(:,t_kernels>0 & t_kernels<0.2,:,:),[],2),4),[1,3,2]),cat(4,temp_plot_tace{:}))
plot_max_error=feval(@(c) permute(std(max(c(:,t_kernels>0 & t_kernels<0.2,:,:),[],2),0,4,'omitmissing')./ sqrt(size(c,4)),[1,3,2]),cat(4,temp_plot_tace{:}))

colors={[1 0 0],[0 0 0],[1 0.5 0.5],[0.5 0.5 0.5]}

% hold on
nexttile
for curr_roi=[ 1 3 ]
    hold on
    ap.errorfill(1:4,plot_max_mean(curr_roi,1:4),plot_max_error(curr_roi,1:4))
end

ylim(0.0001*[0.8 2.2])
xlim([1 4])
xticklabels({'task 1','passive 1','task 2','passive 2'})
legend({'','mPFC','','aPFC'},'Box','off','Location','northeastoutside')


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_H.eps'), ...
    'ContentType','vector');


%%  nose movement in passive
clear all
groups_name={'VA','AV'};
modes_name={'Visual','Auditory'};
U_master = plab.wf.load_master_U;
load('C:\Users\dsong\Documents\MATLAB\Da_Song\DS_scripts_ptereslab\General_information\roi.mat');
Path='D:\Data process\project_cross_model\wf_data\data_package';

surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
t_kernels=[-10:30]/surround_samplerate;
period=find(t_kernels>0&t_kernels<0.2);
surround_window = [-0.5,1];
mousecam_framerate = 30;
face_time = surround_window(1):1/mousecam_framerate:surround_window(2);

all_data=struct;

for curr_group=1:2
    switch curr_group
        case 1
            animals = {'DS007','DS010','AP019','AP021','DS011','AP022'};
        case 2
            animals = {'DS000','DS004','DS015','DS016'};
    end
    for curr_mode=1:2

        temp_data_all=table;
        for curr_animal=1:length(animals)
            preload_vars=who;
            animal=animals{curr_animal};
            data_all=matfile(fullfile(Path,[animal '_all_data.mat']));


            select_id_3{curr_group}=(strcmp([data_all.task_name],'stim_wheel_right_stage1')|...
                strcmp([data_all.task_name],'stim_wheel_right_stage2'))&...
                ~cellfun(@isempty ,data_all.wf_lcr_passive)&...
                ~cellfun(@isempty ,data_all.wf_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_lcr_passive);

            select_id_3{3-curr_group}=(strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume')|...
                strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume'))&...
                ~cellfun(@isempty ,data_all.wf_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.wf_lcr_passive) &...
                ~cellfun(@isempty ,data_all.face_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_lcr_passive);


            switch curr_mode
                case 1
                    temp_face_passive=data_all.face_lcr_passive;
                    temp_wf_passive_0=data_all.wf_lcr_passive;
                case 2
                    temp_face_passive=data_all.face_hml_passive_audio;
                    temp_wf_passive_0=data_all.wf_hml_passive_audio;
            end

            % pupil

            % temp_idd=find(~cellfun(@isempty,temp_face_passive,'UniformOutput',true));
            % pupil_idx=ismember(1:numel(temp_face_passive),temp_idd(cellfun(@(x)  x.validation.pupil, ...
            %     temp_face_passive(temp_idd),'UniformOutput',true)))';


            temp_behavior=data_all.behavior_task;
            temp_p_val=cellfun(@(x) ...
                arrayfun(@(id)  temp_behavior{id}.rxn_l_p(1)<0.05, find(x),'UniformOutput',true),...
                select_id_3,'UniformOutput',false);



            temp_pupil_size=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx)  feval(@(d) cat(3,d{:}), ...
                cellfun(@(c) nanmean(c,1), xx.pupil_data.diameterZ_filt_sav,'uni',false)),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x ),select_id_3,'UniformOutput',false),'UniformOutput',false);



            pupil_data=cat(4,nan(size(temp_pupil_size{1},1), size(temp_pupil_size{1},2),size(temp_pupil_size{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_pupil_size{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_pupil_size{2}(:,:,:,1:min(5,end)));

            pupil_3stage={temp_pupil_size{1}(:,:,:, temp_p_val{1}==0),...
                temp_pupil_size{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_pupil_size{2}(:,:,:, find(temp_p_val{2}==1,2,'last'))};


            %pupil center

            temp_pupil_center=cellfun(@(aa)  feval(@(a) cat(5,a{:}), cellfun(@(xx)  feval(@(d) cat(4,d{:}), ...
                cellfun(@(c) nanmean(c,1), xx.pupil_data.center_filt_sav,'uni',false)),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x),select_id_3,'UniformOutput',false),'UniformOutput',false);



            pupil_center_data=cat(5,nan(size(temp_pupil_center{1},1), size(temp_pupil_center{1},2),...
                size(temp_pupil_center{1},3),size(temp_pupil_center{1},4),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_pupil_center{1}(:,:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_pupil_center{2}(:,:,:,:,1:min(5,end)));


            pupil_center_3stage={temp_pupil_center{1}(:,:,:,:,temp_p_val{1}==0),...
                temp_pupil_center{1}(:,:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_pupil_center{2}(:,:,:,:,find(temp_p_val{2}==1,2,'last'))};



            % nose
            temp_nose_passive=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx) ...
                permute( nanmean(...
                feval(@(d) cat(5,d{:}), cellfun(@(c) nanmean(c(:,:,3,:),1), xx.face_data.nose_filt_sav,'uni',false)) ,1),[2,5,4,1,3]) ,...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x),select_id_3,'UniformOutput',false),'UniformOutput',false);

            nose_data=cat(4,nan(size(temp_nose_passive{1},1), size(temp_nose_passive{1},2),size(temp_nose_passive{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_nose_passive{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_nose_passive{2}(:,:,:,1:5));

            nose_3stage={temp_nose_passive{1}(:,:,:,temp_p_val{1}==0),...
                temp_nose_passive{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_nose_passive{2}(:,:,:,find(temp_p_val{2}==1,2,'last'))};



            % wf_passive
            temp_wf_passive=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx) cat(3,xx.kernels_decoding{:}),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_wf_passive_0(x),select_id_3,'UniformOutput',false),'UniformOutput',false);

            wf_passive_3stage={temp_wf_passive{1}(:,:,:, temp_p_val{1}==0),...
                temp_wf_passive{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_wf_passive{2}(:,:,:, find(temp_p_val{2}==1,2,'last'))};


            wf_passive_data=cat(4,nan(size(temp_wf_passive{1},1), size(temp_wf_passive{1},2),size(temp_wf_passive{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_wf_passive{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_wf_passive{2}(:,:,:,1:5));




            temp_data_all.nose_passive{curr_animal}=nose_data;
            temp_data_all.pupil_size{curr_animal}=pupil_data;
            temp_data_all.pupil_size_single{curr_animal}=pupil_3stage;
            temp_data_all.pupil_all_trace{curr_animal}=temp_pupil_size;

            temp_data_all.wf_passive{curr_animal}=wf_passive_data;
            temp_data_all.wf_passive_single{curr_animal}=wf_passive_3stage;
            temp_data_all.wf_passive_all_trace{curr_animal}=temp_wf_passive;

            temp_data_all.pupil_center{curr_animal}=pupil_center_data;
            temp_data_all.pupil_center_single{curr_animal}=pupil_center_3stage;
            temp_data_all.pupil_center_all_trace{curr_animal}=temp_pupil_center;

            temp_data_all.nose{curr_animal}=nose_data;
            temp_data_all.nose_single{curr_animal}=nose_3stage;
            temp_data_all.nose_all_trace{curr_animal}=temp_nose_passive;


            clearvars('-except',preload_vars{:});
        end

        all_data.([groups_name{curr_group} '_' modes_name{curr_mode}])=temp_data_all;



    end

end


passive_id={[3 1 2],[2 1 3]};
figure('Position',[50 50 1000 800])
mainfig=tiledlayout( 3,3, ...
    'TileSpacing', 'tight', 'Padding', 'none');

plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 1;  % 明确放在主 layout 的第 1 个 tile

% pupil size  trace
colors={[0 0 1;0 0 0],[1 0 0;0 0 0]}
for curr_group=1:2
    for curr_mode=curr_group

        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_size_single{:}));
        temp_pupil_trace=cellfun(@(x) permute(diff(x,1,2),[2,3,4,1]),pupil_data_single,'UniformOutput',false);
        temp_pupil_trace_2=cellfun(@(x) cat(2,x(:,passive_id{curr_mode}(1),:),nanmean(x(:,passive_id{curr_mode}(2:3),:),2)),...
            temp_pupil_trace,'UniformOutput',false);
        % figure('Position',[50 50 400 150]);
        % tiledlayout(1,2)
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_trace_2{curr_state},3),...
                nanstd(temp_pupil_trace_2{curr_state},0,3)./sqrt(size(temp_pupil_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([-0.01 0.03])

            xline(0)
            axis off
        end


    end
end

hold on
line([-0.2 -0.2],[-0.01 0],'Color','k');
line([-0.2 0],[-0.01 -0.01],'Color','k');
text(-0.32, -0.01, '0.01 arb','Rotation',90)
text(-0.2, -0.015, '0.2 s')

% pupil size across days
plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor"  );
plot_fig.Layout.Tile = 2;  % 明确放在主 layout 的第 1 个 tile

colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2
        % pupil
        pupil_data=cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_size{:});
        temp_avg_pupil_data= cat(3,pupil_data(:,:,passive_id{curr_mode}(1),:),nanmean(pupil_data(:,:,passive_id{curr_mode}(2:3),:),3))
        temp_pupil_AUC_mean=feval(@(a)...
            permute(nanmean(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),1),...
            [4,3,2,1]),diff(temp_avg_pupil_data,1,2));

        temp_pupil_AUC_error=feval(@(a)...
            permute(nanstd(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),0,1)./sqrt(size(pupil_data,1)),...
            [4,3,2,1]),diff(temp_avg_pupil_data,1,2));


        nexttile(plot_fig)
        hold on
        ap.errorfill(1:5,temp_pupil_AUC_mean(1:5,:),temp_pupil_AUC_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_pupil_AUC_mean(6:10,:),temp_pupil_AUC_error(6:10,:),[colors_temp{curr_group};0 0 0])

        xlim([1 10])
        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})
        ylim([-0.1 0.4])
        ylabel('Pupil diameter (arb)')

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end


% wf vs pupil size
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)
for curr_group=1:2
    for curr_mode=curr_group

        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_all_trace{:}));

        temp_pupil_AUC_single=cellfun(@(x) feval(@(a)...
            permute(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),...
            [4,3,1,2]),diff(x,1,2)),pupil_data_single,'UniformOutput',false)



        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}))
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_pupil_AUC_single,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Pupil diameter')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(-0.1, 0.7, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end

set(gca,'Color','none')
axis square
drawnow


plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 4;  % 明确放在主 layout 的第 1 个 tile

% pupil saccade  trace
passive_id={[3 1 2],[2 1 3]};
colors={[0 0 1;0 0 0;1 0 0],[1 0 0;0 0 0;1 0 0]}
for curr_group=1:2
    for curr_mode=curr_group




        pupil_center_single=feval(@(a) arrayfun(@(id) cat(5,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center_single{:}));

        % temp_pupil_center_trace=cellfun(@(x) permute( vecnorm((x-x(:,1,:,:,:)),2,3),[2,3,4,5,1]),pupil_center_single,'UniformOutput',false);
        temp_pupil_center_trace=cellfun(@(x) permute( vecnorm(diff(x,1,2),2,3),[2,4,5,3,1]),pupil_center_single,'UniformOutput',false);

        temp_pupil_center_trace_2=cellfun(@(x) cat(2,x(:,passive_id{curr_mode}(1),:),nanmean(x(:,passive_id{curr_mode}(2:3),:),2)),...
            temp_pupil_center_trace,'UniformOutput',false);


        % figure('Position',[50 50 400 150]);
        % tiledlayout(1,2)
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_center_trace_2{curr_state},3),...
                nanstd(temp_pupil_center_trace_2{curr_state},0,3)./sqrt(size(temp_pupil_center_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([0 0.03])
            xline(0)
            axis off
        end


    end
end
hold on
line([-0.2 -0.2],[0 0.01],'Color','k');
line([-0.2 0],[0 0],'Color','k');
text(-0.42, 0.00, '0.01 arb','Rotation',90)
text(-0.2, -0.005, '0.2 s')



plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor" );
plot_fig.Layout.Tile = 5;  % 明确放在主 layout 的第 1 个 tile

% pupil saccade
passive_id={[3 1 2],[2 1 3]};
% fig1=figure;
% tiledlayout(fig1,2,2,"TileIndexing","columnmajor")
colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2

        % pupil center
        pupil_center_data=cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center{:});

        temp_pupil_center=permute(cat(4,vecnorm(diff(pupil_center_data(:,:,:,passive_id{curr_mode}(1),:),1,2),2,3),...
            nanmean(vecnorm(diff(pupil_center_data(:,:,:,passive_id{curr_mode}(2:3),:),1,2),2,3),4)), [2,4,5,1,3]);


        temp_pupil_center_mean=permute(nanmean(trapz(temp_pupil_center(face_time>0&face_time<1,:,:,:),1),4),[3,2,1]);
        temp_pupil_center_error=permute(nanstd(trapz(temp_pupil_center(face_time>0&face_time<1,:,:,:),1),0,4)./...
            sqrt(size(temp_pupil_center,4)),[3,2,1]);

        nexttile(plot_fig)
        ap.errorfill(1:5,temp_pupil_center_mean(1:5,:),temp_pupil_center_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_pupil_center_mean(6:10,:),temp_pupil_center_error(6:10,:),[colors_temp{curr_group};0 0 0])
        ylim([0 1])
        xlim([1 10])
        ylabel('Pupil saccade (arb)')

        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end

% wf vs pupil saccade
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 6;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)

for curr_group=1:2
    for curr_mode=curr_group

        pupil_center_single=feval(@(a) arrayfun(@(id) cat(5,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center_all_trace{:}));

        temp_pupil_center=cellfun(@(x) permute(vecnorm(diff(x,1,2),2,3),[5,4,2,1,3]),...
            pupil_center_single,'UniformOutput',false);

        temp_pupil_center_mean=cellfun(@(x) trapz(x(:,:,face_time>0&face_time<1),3),temp_pupil_center,'UniformOutput',false);

        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}));
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_pupil_center_mean,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Pupil saccade')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(0, 1, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end

set(gca,'Color','none')
axis square
drawnow



% nose move  trace

plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 7;  % 明确放在主 layout 的第 1 个 tile
colors={[0 0 1;0 0 0],[1 0 0;0 0 0]}
for curr_group=1:2
    for curr_mode=curr_group
        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose_single{:}));
        temp_pupil_trace=cellfun(@(x) permute(vecnorm(diff(x,1,1),2,3),[2,1,4,3]),pupil_data_single,'UniformOutput',false);
        temp_pupil_trace_2=cellfun(@(x) cat(1,x(passive_id{curr_mode}(1),:,:),nanmean(x(passive_id{curr_mode}(2:3),:,:),1)),...
            temp_pupil_trace,'UniformOutput',false);
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_trace_2{curr_state},3)',...
                nanstd(temp_pupil_trace_2{curr_state},0,3)'./sqrt(size(temp_pupil_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([-0.5 1])
            xline(0)
            axis off
        end
    end
end

hold on
line([-0.2 -0.2],[-0.5 0],'Color','k');
line([-0.2 0],[-0.5 -0.5],'Color','k');
text(-0.42,  -0.4, '0.5 arb' ,'Rotation',90)
text(-0.1, -0.6, '0.2 s')

% nose move across days
passive_id={[3 1 2],[2 1 3]};

plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor" );
plot_fig.Layout.Tile = 8;  % 明确放在主 layout 的第 1 个 tile

colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2
        % pupil
        nose_data=cat(5, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose{:});
        temp_avg_nose_data= cat(2,vecnorm(diff(nose_data(:,passive_id{curr_mode}(1),:,:,:),1,1),2,3),...
            nanmean(vecnorm(diff(nose_data(:,passive_id{curr_mode}(2:3),:,:,:),1,1),2,3),2));
        temp_nose_AUC_mean= permute(nanmean(trapz(1:sum(face_time>0&face_time<1),temp_avg_nose_data(face_time>0&face_time<1,:,:,:,:),1),5),...
            [4,2,3,1]);

        temp_nose_AUC_error=...
            permute(nanstd(trapz(1:sum(face_time>0&face_time<1),temp_avg_nose_data(face_time>0&face_time<1,:,:,:,:),1),0,5)./sqrt(size(nose_data,5)),...
            [4,2,3,1]);


        nexttile(plot_fig)
        hold on
        ap.errorfill(1:5,temp_nose_AUC_mean(1:5,:),temp_nose_AUC_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_nose_AUC_mean(6:10,:),temp_nose_AUC_error(6:10,:),[colors_temp{curr_group};0 0 0])

        xlim([1 10])
        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})
        ylim([0 15])
        ylabel('Nose move (arb)')

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end

% nose move vs wf
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 9;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)

for curr_group=1:2
    for curr_mode=curr_group

        nose_move=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose_all_trace{:}));

        nose_move_trace=cellfun(@(x) permute(vecnorm(diff(x,1,1),2,3),[4,2,1,3]),...
            nose_move,'UniformOutput',false);

        nose_move_trace_mean=cellfun(@(x) trapz(x(:,:,face_time>0&face_time<1),3),nose_move_trace,'UniformOutput',false);

        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}));
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            nose_move_trace_mean,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Nose move')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(0, 20, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end
set(gca,'Color','none')

axis square
drawnow


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_I.eps'), ...
    'ContentType','vector');




%%    EDF J    mPFC aPFC ephys 
clear all

Path='D:\Data process\project_cross_model\wf_data\data_package';
U_master = plab.wf.load_master_U;

surround_samplerate = 35;
surround_window_passive = [-0.5,1];
surround_window_task = [-0.2,1];
t_kernels=1/surround_samplerate*[-10:30];
period_kernels=find(t_kernels>0&t_kernels<0.2);



animals={'DS029','DS030','DS031'};
temp_kernels=cell(length(animals),1);

for curr_animal=1:length(animals)
    animal=animals{curr_animal};
    temp_path=matfile(fullfile(Path,[animal '_all_data.mat']));
    temp_name_idx=cellfun(@(x) strcmp(x,'stim_wheel_right_stage2_mixed_VA_earphone'),temp_path.task_name,'uni',true );
    temp_task=temp_path.wf_task;
     temp_kernels{curr_animal}= temp_task{~cellfun(@isempty ,temp_task,'UniformOutput',true)&temp_name_idx}.stim_kernels;
end

wf_task=...
feval(@(aa,bb)    cat(4,aa,bb),...
feval(@(c)  cat(3,c{:}) ,arrayfun(@(id)  temp_kernels{id}{1}   , 1:length(animals),'UniformOutput',false)),...
feval(@(c)  cat(3,c{:}) ,arrayfun(@(id)  temp_kernels{id}{2}   , 1:length(animals),'UniformOutput',false)));

temp_image=plab.wf.svd2px(U_master(:,:,1:size(wf_task,1)),wf_task);


image_max=permute(nanmean(max(temp_image(:,:,period_kernels,:,:),[],3),4),[1,2,5,4,3]);

figure('Position',[50 50 300 600])
Color={'B','R'}
probe_colors={[0.5 1 0.5],[0.2 0.5 0.2],[0 0.8 0]}
h=struct
for curr_mod=1:2
    ax=nexttile
    imagesc(ax,image_max(:,:,curr_mod))
    axis image off;
    clim(0.0003 .* [0, 1]);
    colormap(ax, ap.colormap(['W' Color{curr_mod}] ));
    hold on
    for curr_animal=1:length(animals)

        ap.wf_draw('probe_insertion',animals{curr_animal});
        h.(['l' num2str(curr_animal)]) = findobj(gca,'Type','line');
        for i = 1:length( h.(['l' num2str(curr_animal)]))
            if strcmp(h.(['l' num2str(curr_animal)])(i).Marker,'.')&h.(['l' num2str(curr_animal)])(i).MarkerSize==20
                h.(['l' num2str(curr_animal)])(i).Color = probe_colors{curr_animal};
                h.(['l' num2str(curr_animal)])(i).MarkerSize = 10;
            end
        end
    end
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_Ja0.eps']), ...
    'ContentType','vector'); 



Path='D:\Data process\project_cross_model\wf_data\data_package';
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');
period=t_bins>0&t_bins<0.2;

animals={'DS029','DS030','DS031'};

% mPFC_idx={[2,4],[2,4,5],[3,4,5]}

responsive_visual_all=cell(2,1);
responsive_audio_all=cell(2,1);
responsive_or_all=cell(2,1);
visual_passive_all=cell(2,1);
audio_passive_all=cell(2,1);
visual_task_all=cell(2,1);
audio_task_all=cell(2,1);
ephys_data_all=cell(2,1);
for curr_group=1:2
    switch curr_group
        case 1
            PFC_idx={[1 ],[3],[1 2]};
        case 2
            PFC_idx={[2,4],[2,4,5],[3,4,5]};
    end


    visual_task=cell(length(animals),1);
    audio_task=cell(length(animals),1);
    visual_passive=cell(length(animals),1);
    audio_passive=cell(length(animals),1);
    responsive_idx=cell(length(animals),1);
    responsive_visual=cell(length(animals),1);
    responsive_audio=cell(length(animals),1);
    temp_ephys=cell(length(animals),1);
    for curr_animal=1:length(animals)
        animal=animals{curr_animal};
        temp_path=matfile(fullfile(Path,[animal '_all_data.mat']));
        % temp_path.data_all_index


        temp_audio_passive=temp_path.ephys_hml_passive_audio_earphone;
        temp_visual_passive=temp_path.ephys_lcr_passive;
        temp_task=temp_path.ephys_task;
        temp_behavior=temp_path.behavior_task;
        ephys_idx=~cellfun(@isempty ,temp_visual_passive );

        ephys_audio=temp_audio_passive(ephys_idx);
        ephys_visual=temp_visual_passive(ephys_idx);
        ephys_task=temp_task(ephys_idx);
        behavior=temp_behavior(ephys_idx);




        for curr_day=1:length(PFC_idx{curr_animal})
            temp_day=PFC_idx{curr_animal}(curr_day);


            response_idx=ephys_visual{temp_day}.response_p{3}>0.95|ephys_audio{temp_day}.response_p{2}>0.95;
            depth_idx=ephys_task{temp_day}.depth>1840;
            responsive_idx{curr_animal}{curr_day}=response_idx(find(depth_idx));
            responsive_visual{curr_animal}{curr_day}=feval(@(a) a((find(depth_idx))), ephys_visual{temp_day}.response_p{3}>0.95);
            responsive_audio{curr_animal}{curr_day}=feval(@(a) a((find(depth_idx))), ephys_audio{temp_day}.response_p{2}>0.95);

            visual_task{curr_animal}{curr_day}= ephys_task{temp_day}.psth(depth_idx,:,1);
            audio_task{curr_animal}{curr_day}= ephys_task{temp_day}.psth(depth_idx,:,2);
            visual_passive{curr_animal}{curr_day}=ephys_visual{temp_day}.psth(depth_idx,:,3);
            audio_passive{curr_animal}{curr_day}=ephys_audio{temp_day}.psth(depth_idx,:,2);

            temp_ephys{curr_animal}{curr_day}= cat(3,ephys_task{temp_day}.psth(depth_idx,:,:),...
                ephys_visual{temp_day}.psth(depth_idx,:,:),...
                 ephys_audio{temp_day}.psth(depth_idx,:,:));

        end


    end


    responsive_visual_all{curr_group}=cat(2,responsive_visual{:});
    responsive_audio_all{curr_group}=cat(2,responsive_audio{:});
    responsive_or_all{curr_group}=cat(2,responsive_idx{:});
    visual_passive_all{curr_group}=cat(2,visual_passive{:});
    audio_passive_all{curr_group}=cat(2,audio_passive{:});
    visual_task_all{curr_group}=cat(2,visual_task{:});
    audio_task_all{curr_group}=cat(2,audio_task{:});
    ephys_data_all{curr_group}=cat(2,temp_ephys{:});

end

temp_or_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_or_all,'uni',false);
temp_v_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_visual_all,'uni',false);
temp_a_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_audio_all,'uni',false);
temp_v_frac=cellfun(@(x) cellfun(@(a) sum(a)/length(a) ,x,'UniformOutput',true) ,responsive_visual_all,'UniformOutput',false   );
temp_a_frac=cellfun(@(x) cellfun(@(a) sum(a)/length(a) ,x,'UniformOutput',true) ,responsive_audio_all,'UniformOutput',false   );
temp_frac=cellfun(@(x,y)    {x,y}, temp_v_frac,temp_a_frac,'UniformOutput',false   );


idx_A_only = cellfun(@(x,y)  find(x&~y),temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_AB     = cellfun(@(x,y) find (x&y),temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_B_only = cellfun(@(x,y)  find(~x&y),temp_v_idx,temp_a_idx,'UniformOutput',false);

idx_order=cellfun(@(x,y)  [find(x&~y); find(x&y); find(~x&y)],temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_orders=cellfun(@(x,y)  {find(x&~y); find(x&y); find(~x&y)},temp_v_idx,temp_a_idx,'UniformOutput',false);

temp_v_p=cellfun(@(x)     cat(1,x{:}), visual_passive_all,'UniformOutput',false);
temp_a_p=cellfun(@(x)     cat(1,x{:}), audio_passive_all,'UniformOutput',false);
temp_v_t=cellfun(@(x)     cat(1,x{:}), visual_task_all,'UniformOutput',false);
temp_a_t=cellfun(@(x)     cat(1,x{:}), audio_task_all,'UniformOutput',false);

temp_all=cellfun(@(x)     cat(1,x{:}), ephys_data_all,'UniformOutput',false);



temp_v_p_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    visual_passive_all,responsive_or_all,'UniformOutput',false);
temp_a_p_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    audio_passive_all,responsive_or_all,'UniformOutput',false);

temp_v_t_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    visual_task_all,responsive_or_all,'UniformOutput',false);
temp_a_t_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    audio_task_all,responsive_or_all,'UniformOutput',false);



for curr_group=1:2
    figure
    tiledlayout(2,11,"TileIndexing","columnmajor")
    for curr_stat=1:11
        nexttile
        imagesc(t_bins,[],temp_all{curr_group}(idx_order{curr_group},:,curr_stat) )
        clim([0 3])
        xlim([-0.2 0.5])
        colormap(ap.colormap(['WB']))
        yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
        yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
        nexttile
        hold on
        for curr_type=1:3
            ap.errorfill(t_bins,nanmean(temp_all{curr_group}(idx_orders{curr_group}{curr_type},:,curr_stat),1),...
                std(temp_all{curr_group}(idx_orders{curr_group}{curr_type},:,curr_stat),0,1)./100)
            ylim([0 5])
        end
    end
end

for curr_group=1:2


    figure;
    mainfig=tiledlayout(1,3,'TileSpacing','tight')

    % A1=nexttile
    % imagesc(t_bins,[],temp_v_t{curr_group}(idx_order{curr_group},:) )
    % clim([0 3])
    % xlim([-0.2 0.5])
    % colormap(A1,ap.colormap(['WB']))
    % yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    % yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    % xline(0)
    % axis off
    % A2=nexttile
    % imagesc(t_bins,[],temp_a_t{curr_group}(idx_order{curr_group},:) )
    % clim([0 3])
    % xlim([-0.2 0.5])
    % colormap(A2,ap.colormap(['WR']))
    % yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    % yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    % xline(0)
    % axis off

    A1=nexttile
    imagesc(t_bins,[],temp_v_p{curr_group}(idx_order{curr_group},:) )
    clim([0 3])
    xlim([-0.2 0.5])
    colormap(A1,ap.colormap(['WB']))
    yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    xline(0)

    axis off
    A2=nexttile
    imagesc(t_bins,[],temp_a_p{curr_group}(idx_order{curr_group},:) )
    clim([0 3])
    xlim([-0.2 0.5])
    colormap(A2,ap.colormap(['WR']))
    yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    xline(0)
    axis off


    plot_fig=tiledlayout(mainfig,2 ,1, ...
        'TileSpacing', 'tight');
    plot_fig.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile
    a4=nexttile(plot_fig,1)

    hold on
    ap.errorfill(t_bins,nanmean(temp_v_p_mean{curr_group},1),nanstd(temp_v_p_mean{curr_group},0,1)./sqrt(size(temp_v_p_mean{curr_group},1)),[0 0 1])
    ap.errorfill(t_bins,nanmean(temp_a_p_mean{curr_group},1),nanstd(temp_a_p_mean{curr_group},0,1)./sqrt(size(temp_a_p_mean{curr_group},1)),[1 0 0])
    ylim([-0.1 2.5])
    xlim([-0.2 0.5])
    xline(0)
    axis off

    a4=nexttile(plot_fig,2)
    ds.make_bar_plot(temp_frac{curr_group},'ColorCell',{[0 0 1],[1 0 0]})
    ylabel('Fraction')
    ylim([0 0.3])
    yticks([0 0.3])
    xticks([])
    set(gca,'Color','none')


ds.shuffle_test(temp_frac{curr_group}{1},temp_frac{curr_group}{2},1)


exportgraphics(gcf, fullfile(plab.locations.server_path,...
    ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_Ja' num2str(curr_group)   '.eps']), ...
    'ContentType','vector'); 

end

response_sort_both=...
    cellfun(@(a,b,c)     cellfun(@(x,y,z)   [x(z) y(z)] ,a,b,c,'uni',false)',...
    responsive_visual_all,responsive_audio_all,...
    responsive_or_all,'UniformOutput',false);

response_sort_all_both=cellfun(@(x)   cat(1,x{:}),response_sort_both,'UniformOutput',false);
response_both_fration_both = cellfun(@(a)  cellfun(@(x) sum(sum(x,2)==2)./size(x,1),a,'UniformOutput',true) ,...
    response_sort_both,'UniformOutput',false);

n_shuff=1000;
response_both_shuff_both= cellfun(@(x)  arrayfun(@(shuff) ap.shake(x,1),1:n_shuff,...
    'UniformOutput',false),response_sort_all_both,'UniformOutput',false);


response_both_fration_shuff_both=cellfun(@(a)  cellfun(@(x) sum( sum(x,2)==2)/ length(x),a,'UniformOutput',true),...
    response_both_shuff_both,'UniformOutput',false);
line_mean_both=cellfun(@(a)  (prctile(a,95)+ ...
    prctile(a,5))/2,response_both_fration_shuff_both,'UniformOutput',false);
line_error_both=cellfun(@(a) (prctile(a,95)-...
    prctile(a,5))/2,response_both_fration_shuff_both,'UniformOutput',false);







colors = {[0 0 1],[ 1.0, 0.647, 0.0],[1 0 0]};
title_names={'mPFC','aPFC'}

scales_all=[0 20 80];
scales = {[0 20 80],[0 8 20],[0 4 10]};
colors = {[0 0 1],[ 1.0, 0.647, 0.0],[1 0 0]};
for curr_group=1:2
    figure('Position',[50 50 200 200])
    ax1 = axes('Position',[0.2 0.2 0.5 0.5]); % 左边大一些
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    xlim(scales_all(1:2))
    ylim(scales_all(1:2))
    xticks(scales_all(1:2))
    yticks(scales_all(1:2))
    set(gca,'Color','none')
    xlabel({'visual response' ;'(\DeltaFR/FR_{0})'})
    ylabel({'auditory response';' (\DeltaFR/FR_{0})'})

    ax2 = axes('Position',[0.75 0.2 0.2 0.5]);
    hold on

    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    xlim(scales_all(2:3))
    ylim(scales_all(1:2))
    xticks(scales_all(3))
    set(gca,'YTick',[])
    set(gca,'YColor','none')
    set(gca,'Color','none')

    ax3 = axes('Position',[0.2 0.75 0.5 0.2]);
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    ylim(scales_all(2:3))
    xlim(scales_all(1:2))
    set(gca,'Color','none')
    yticks(scales_all(3))
    set(gca,'XTick',[])
    set(gca,'XColor','none')

    ax4 = axes('Position',[0.75 0.75 0.2 0.2]);
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    ylim(scales_all(2:3))
    xlim(scales_all(2:3))
    axis off
    title(title_names{curr_group},'FontWeight','normal')

    annotation('line',[0.74 0.76],[0.19 0.21],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.69 0.71],[0.19 0.21],'Color','k','LineWidth',1) % 右斜杠
    annotation('line',[0.19 0.21],[0.74 0.76],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.19 0.21],[0.69 0.71],'Color','k','LineWidth',1) % 右斜杠

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_Jf' num2str(curr_group)   '.eps']), ...
    'ContentType','vector'); 

end



figure('Position',[50 50 100 200])
colors={[0 0.5 0],[0.5 0 0.2]};
maker_size=20;

hold on
ds.make_bar_plot(response_both_fration_both,'ColorCell',colors,'BarAlpha',0.2,'DotSize',maker_size)
arrayfun(@(curr_group)  ap.errorfill(curr_group-0.4:0.8:curr_group+0.4,[line_mean_both{curr_group} line_mean_both{curr_group}],...
    [line_error_both{curr_group} line_error_both{curr_group}],[0.5 0.5 0.5],1),1:2);
ylim([0 0.5])
yticks([0 0.5])
xticks([1 2])
xticklabels({'mPFC','aPFC'})
ylabel(['(V∩A)/(V∪A)' ])
set(gca,'color','none')

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_Jg' num2str(curr_group)   '.eps']), ...
    'ContentType','vector'); 
