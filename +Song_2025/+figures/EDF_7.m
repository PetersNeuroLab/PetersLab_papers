%% EDF 7   viarable task  contrast volume
main_preload_vars = who;


tempdata{1}= load(fullfile(Path,'data\revision\visual_task_variable.mat'));
tempdata{2}=load(fullfile(Path,'data\revision\audio_task_variable.mat'));


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
    temp_image_max=permute(nanmean(max(temp_images(:,:,kernels_period,:,:),[],3),5),[1,2,4,3,5]);
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
    temp_plot_tace_max=permute(max(temp_plot_tace(:,kernels_period,:,:),[],2),[3,1,4,2]);

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
    temp_plot_tace_max=permute(max(temp_plot_tace(:,kernels_period,:,:),[],2),[3,1,4,2]);

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
 clearvars('-except',main_preload_vars{:});
