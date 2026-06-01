%% Fig s15 b,c,d,e behavioral results in mixed tasks
main_preload_vars = who;

load(fullfile(Path,'data','behavior'));

perform_mix=cellfun(@(x) structfun(@(a) [nanmean(a(end-5 :end-3,:),1); nanmean(a(end-2 :end,:),1)]' ,...
    x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false);
perform_mix_mean=cellfun(@(x) structfun(@(a) mean([nanmean(a(end-5 :end-3,:),1); nanmean(a(end-2 :end,:),1)]',1,'omitmissing') ,...
    x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false);

perform_mix_error=cellfun(@(x) structfun(@(a) std([nanmean(a(end-5 :end-3,:),1); nanmean(a(end-2 :end,:),1)]',0,1,'omitmissing')./sqrt(size(a,2)) ,...
    x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false);


% perform_mix_mean{curr_group}.reaction_time
barColors =  [[84 130 53]./255; [112  48 160 ]./255]; % 浅蓝、浅红
scatterColors = [[84 130 53]./255; [112  48 160 ]./255]; % 深蓝、深红

figure('Position', [50 50 500 150]);
t1 = tiledlayout(1,4, 'TileSpacing', 'loose', 'Padding', 'loose');
nexttile()
hold on
for curr_group=1:2
    for    curr_task=[0 2]
        switch curr_task
            case 0
                ii=1
            case 2
                ii=2
        end

        errorbar(curr_task+curr_group, perform_mix_mean{curr_group}.reaction_time(ii),...
            perform_mix_error{curr_group}.reaction_time(ii),...
            'o','LineStyle', 'none',...
            'CapSize', 0,...
            'MarkerEdgeColor',scatterColors(curr_group,:) , ...
            'MarkerFaceColor',barColors(curr_group,:) , ...
            'Color', barColors(curr_group,:),...
            'LineWidth',1.5,'MarkerSize',4.5)
    end

end

p_val_corss_group= arrayfun(@(id) ds.shuffle_test( perform_mix{1}.reaction_time(:,id),perform_mix{2}.reaction_time(:,id),0,1),1:2,'UniformOutput',true)
p_val_within_group= arrayfun(@(id) ds.shuffle_test( perform_mix{id}.reaction_time(:,1),perform_mix{id}.reaction_time(:,2),1,1),1:2,'UniformOutput',true)

xlim([0 5])
ylabel('RT (s)')
ylim([0 0.5])
yticks([0 0.5])
xticks([1.5 3.5])
xticklabels({'mixed V','mixed A'})
% set(gca, 'YScale', 'log', 'Color', 'none');
set(gca, 'Color', 'none');


for curr_state=1:2
    if p_val_within_group(curr_state) < 0.05
        stars = repmat('*',1,sum(p_val_within_group(curr_state)<[0.05 0.01 0.001]));
        y_sig = 0.2+0.1*curr_state;
        plot([curr_state curr_state+2], [1 1]*y_sig, 'k-');
        text(curr_state+1, y_sig+0.02*curr_state, stars, 'HorizontalAlignment','center');
    end
end
nexttile
hold on
for curr_group=1:2
    errorbar([0 2]+curr_group, perform_mix_mean{curr_group}.performance,...
        perform_mix_error{curr_group}.performance,...
        'o','LineStyle', 'none',...
        'CapSize', 0,...
        'MarkerEdgeColor', scatterColors(curr_group,:), ...
        'MarkerFaceColor', scatterColors(curr_group,:), ...
        'Color', scatterColors(curr_group,:),...
        'LineWidth',1.5,'MarkerSize',4.5)
end
p_val_corss_group= arrayfun(@(id) ds.shuffle_test( perform_mix{1}.performance(:,id),perform_mix{2}.performance(:,id),0,1),1:2,'UniformOutput',true)
p_val_within_group= arrayfun(@(id) ds.shuffle_test( perform_mix{id}.performance(:,1),perform_mix{id}.performance(:,2),0,1),1:2,'UniformOutput',true)


ylabel('performance')
xlim([0 5])
ylim([0 0.8])
xticks([1.5 3.5])
xticklabels({'mixed V','mixed A'})
yticks([0 0.8])
set(gca, 'Color', 'none');


nexttile
hold on
for curr_group=1:2
    hold on
    errorbar(curr_group,perform_mix_mean{curr_group}.itimove(2),perform_mix_error{curr_group}.itimove(2),...
        'o','Color',scatterColors(curr_group,:),'LineWidth',1.5,...
        'MarkerFaceColor',scatterColors(curr_group,:),'MarkerSize',4.5,'CapSize',0);
    % scatter(curr_group,itimove_temp_mix_peak_mean{curr_group},'MarkerFaceColor',scatterColors(curr_group,:),'MarkerEdgeColor','none');
end
ylabel('relative move')
xlim([0.5 2.5])
ylim([0 5])
yticks([0 5])

xticks([1.5])
xticklabels({'mixed'})
set(gca, 'Color', 'none');
p_val_corss_group= arrayfun(@(id) ds.shuffle_test( perform_mix{1}.itimove(:,id),perform_mix{2}.itimove(:,id),0,1),1:2,'UniformOutput',true)


nexttile
hold on
for curr_group=1:2
    errorbar([0 2]+curr_group, perform_mix_mean{curr_group}.velocity,...
        perform_mix_error{curr_group}.velocity,...
        'o','LineStyle', 'none',...
        'CapSize', 0,...
        'MarkerEdgeColor', scatterColors(curr_group,:), ...
        'MarkerFaceColor', scatterColors(curr_group,:), ...
        'Color', scatterColors(curr_group,:),...
        'LineWidth',1.5,'MarkerSize',4.5)
end

p_val_corss_group= arrayfun(@(id) ds.shuffle_test( perform_mix{1}.velocity(:,id),perform_mix{2}.velocity(:,id),0,1),1:2,'UniformOutput',true)

p_val_within_group= arrayfun(@(id) ds.shuffle_test( perform_mix{id}.velocity(:,1),perform_mix{id}.velocity(:,2),0,1),1:2,'UniformOutput',true)

for curr_state=1:2
    if p_val_corss_group(curr_state) < 0.05
        stars = repmat('*',1,sum(p_val_corss_group(curr_state)<[0.05 0.01 0.001]));
        y_sig = 3700;
        plot(2*curr_state-1:2*curr_state, [1 1]*y_sig, 'k-');
        text(2*curr_state-0.5, y_sig+100, stars, 'HorizontalAlignment','center');
    end
end

for curr_state=1:2
    if p_val_within_group(curr_state) < 0.05
        stars = repmat('*',1,sum(p_val_corss_group(curr_state)<[0.05 0.01 0.001]));
        y_sig = 3800+curr_state*100;
        plot([curr_state curr_state+2], [1 1]*y_sig, 'k-');
        text(curr_state+1, y_sig+100, stars, 'HorizontalAlignment','center');
    end
end


ylabel('velocity')
xlim([0 5])
ylim([0 4200])
xticks([1.5 3.5])
xticklabels({'mixed V','mixed A'})
yticks([0 3800])
yticklabels({'0','max'})
set(gca, 'Color', 'none');


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s4b.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});

%% fig s15 f

main_preload_vars = who;

load(fullfile(Path,'data','wf_task_kernels'));

tem_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,21:26,:)),  wf_task_kernels_across_day,'UniformOutput',false);
image_max=cellfun(@(x)    cat(3,  nanmean(max(x(:,:,kernels_period,[1: 3],:),[],3),[4 5]),...
    nanmean(max(x(:,:,kernels_period,[4: 6],:),[],3),[4 5])),tem_image,'UniformOutput',false  );

image_max_peak=cellfun(@(x)  ds.make_each_roi( permute( cat(3,  nanmean(max(x(:,:,kernels_period,[1: 3],:),[],3),4),...
    nanmean(max(x(:,:,kernels_period,[4 :6],:),[],3),4)),[1,2,3,5,4]),2,roi1),tem_image,'UniformOutput',false  );

image_max_peak_mean=cellfun(@(x)  nanmean(ds.make_each_roi( permute( cat(3,  nanmean(max(x(:,:,kernels_period,[1: 3],:),[],3),4),...
    nanmean(max(x(:,:,kernels_period,[4: 6],:),[],3),4)),[1,2,3,5,4]),2,roi1),3),tem_image,'UniformOutput',false  );

image_max_peak_error=cellfun(@(x)  std(ds.make_each_roi( permute( cat(3,  nanmean(max(x(:,:,kernels_period,[1: 3],:),[],3),4),...
    nanmean(max(x(:,:,kernels_period,[4: 6],:),[],3),4)),[1,2,3,5,4]),2,roi1),0,3,'omitmissing')./sqrt(size(x,5)),tem_image,'UniformOutput',false  );


figure('Position',[50 50 400 300])
imagelayout = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');

colors = {'G','P'};
titles = {'mixed V','mixed A'};
axs = gobjects(2,2);

% —— 左侧 2×2 图像 —— %
for curr_group = 1:2
    for curr_stage = 1:2
        ax = nexttile(imagelayout, 3*curr_group-3+curr_stage);
        axs(curr_group,curr_stage) = ax;
        imagesc(ax, image_max{curr_group}(:,:,curr_stage));
        axis(ax,'image','off');
        clim(ax, 0.0003*[0 1]);
        ap.wf_draw('ccf',[0.5 0.5 0.5]);
        colormap(ax, ap.colormap(['W' colors{curr_group}]));
        if curr_group==1
            title(ax, titles{curr_stage},'FontSize',10,'FontWeight','normal');
        end
    end
end

% —— 每行横向 colorbar（底下，宽=子图宽的1/3） —— %
drawnow
gap = 0.01; cbw = 0.02; shrink = 1/3;
for r = 1:2
    p = axs(r,2).Position;                      % 右侧子图的位置
    cb = colorbar(axs(r,2),'southoutside');     % 关联右侧子图
    cb.Units = 'normalized';
    w = p(3)*shrink; x0 = p(1)+(p(3)-w)/2; y0 = p(2)+0.01;
    cb.Position = [x0 y0 w cbw];
    cb.Ticks = [cb.Limits(1), cb.Limits(2)]; % 只保留最小值和最大值

    cb.Label.String='\Delta F/F_{0}';
    % ticks = cb.Ticks; cb.TickLabels = arrayfun(@(x)sprintf('%.2f',x),ticks,'uni',false);
end

% —— 右列两幅误差散点图（tile 3 / tile 6） —— %
scale_plot = 0.0005;
categories = {'mixed V','mixed A'};
x = 1:numel(categories);
colors1 = [[84 130 53]; [112 48 160]]/255;

for curr_area = [1 3]
    ax1 = nexttile(imagelayout, 3*ceil(curr_area/2));  % 1->3, 3->6
    hold(ax1,'on');

    arrayfun(@(id) ds.shuffle_test(  image_max_peak{1}(curr_area,id,:),image_max_peak{2}(curr_area,id,:),0,2),1:2,'UniformOutput',true)


    data   = [image_max_peak_mean{1}(curr_area,:);  image_max_peak_mean{2}(curr_area,:)]';
    errors = [image_max_peak_error{1}(curr_area,:); image_max_peak_error{2}(curr_area,:)]';
    offsets = linspace(-0.5,0.5,size(data,2));

    for curr_group = 1:size(data,2)
        x_offset = x*2 + offsets(curr_group) - 0.5;
        errorbar(ax1, x_offset, data(:,curr_group), errors(:,curr_group), ...
            'o','LineWidth',1.5,'MarkerSize',4, ...
            'Color',colors1(curr_group,:), ...
            'MarkerFaceColor',colors1(curr_group,:), ...
            'LineStyle','none','CapSize',0);
    end

    set(ax1,'XTick',x*2-0.5,'XTickLabel',categories,'XLim',[0 5], ...
        'YLim',[0 scale_plot],'Box','off','Color','none');
    ylabel(ax1,'\Delta F/F_{0}');
    title(ax1, roi1(curr_area).name,'FontWeight','normal');

    % inset（左上角 ROI mask）
    mainPos = ax1.Position;
    insetAx = axes('Position', [mainPos(1)+0.05*mainPos(3), mainPos(2)+0.65*mainPos(4), 0.3*mainPos(3), 0.3*mainPos(4)]);
    imagesc(insetAx, roi1(curr_area).data.mask); ap.wf_draw('ccf',[0.5 0.5 0.5]);
    axis(insetAx,'image','off'); ylim(insetAx,[0 200]); xlim(insetAx,[20 220]);
    clim(insetAx,[0 1]); colormap(insetAx, ap.colormap('WK')); uistack(insetAx,'bottom');
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s4c.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});


%% fig s15 g

main_preload_vars = who;

load(fullfile(Path,'data','wf_passive_kernels'));

tem_image=cellfun(@(x) cellfun(@(a) plab.wf.svd2px(U_master(:,:,1:size(a,1)),a(:,:,:,24:26,:)),...
    x,'UniformOutput',false),wf_passive_kernels_across_day,'UniformOutput',false);

image_max=cellfun(@(x)  cellfun(@(a) permute( nanmean(max(a(:,:,kernels_period,:,:,:),[],3),[5 6]),[1,2,4,3]),...
    x,'UniformOutput',false), tem_image,'UniformOutput',false  );

image_max_peak=cellfun(@(x) cellfun(@(a) ds.make_each_roi( permute( nanmean(max(a(:,:,kernels_period,:,:,:),[],3),5),...
    [1,2,4,6,5,3]),2,roi1),x,'UniformOutput',false),tem_image,'UniformOutput',false  );

image_max_peak_mean=cellfun(@(x) cellfun(@(a)  nanmean( ds.make_each_roi( permute( nanmean(max(a(:,:,kernels_period,:,:,:),[],3),5),...
    [1,2,4,6,5,3]),2,roi1),3),x,'UniformOutput',false),tem_image,'UniformOutput',false  );

image_max_peak_error=cellfun(@(x) cellfun(@(a) std( ds.make_each_roi( permute( nanmean(max(a(:,:,kernels_period,:,:,:),[],3),5),...
    [1,2,4,6,5,3]),2,roi1),0,3,'omitmissing')./sqrt(size(a,6)),x,'UniformOutput',false),tem_image,'UniformOutput',false  );




figure('Position',[50 50 400 300])
imagelayout = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');

% —— 左侧 2×2 图像 —— %
axs = gobjects(2,2);
titles = {'V passive','A passive'};
used_stim = [3 2];
CLim_im = 0.0003*[0 1];
colors_lr = {'G','P'};   % 仅用于左侧图像的 colormap 码

for curr_group = 1:2
    for curr_stage = 1:2
        ax = nexttile(imagelayout, 3*curr_group-3+curr_stage);   % (r,c) -> 线性索引
        axs(curr_group,curr_stage) = ax;
        imagesc(ax, image_max{curr_stage}{curr_group}(:,:,used_stim(curr_stage)));
        axis(ax,'image','off');
        clim(ax, CLim_im);
        ap.wf_draw('ccf',[0.5 0.5 0.5]);
        colormap(ax, ap.colormap(['W' colors_lr{curr_group}]));
        if curr_group==1, title(ax, titles{curr_stage},'FontSize',10,'FontWeight','normal'); end
    end
end

% —— 每行一个横向 colorbar（位于该行右图底下；宽=子图宽的1/3） —— %
drawnow
gap = 0.01; cbh = 0.02; shrink = 1/3;
for r = 1:2
    p  = axs(r,2).Position;                         % 右侧子图的位置
    cb = colorbar(axs(r,2),'southoutside');         % 关联右侧子图
    cb.Units = 'normalized';
    w = p(3)*shrink; x0 = p(1)+(p(3)-w)/2; y0 = p(2)-gap-cbh;
    cb.Position = [x0 y0 w cbh];
    % 如需科学计数法刻度，取消注释下一行：
    % cb.TickLabels = arrayfun(@(x) sprintf('%.2e',x), cb.Ticks, 'uni', false);
end

% —— 右列两幅误差散点图（tile 3 / 4 属于行1；tile 7 / 8 属于行2） —— %
colors = {  [84 130 53]./255 ; [112 48 160]./255 };


scale_plot = 0.0003;                % 对应你的 ylim(0.0002*[0 1])
cats = {{'R','C','L'};{'8k','4k','12k'}};
p_val_corss_group=cell(2,1)
% 只放中间一个刻度
p_val=cell(3,1);
for curr_roi=[1 3]
    switch curr_roi
        case 1
            a2 = nexttile(imagelayout, 3);


        case 3
            a2 = nexttile(imagelayout, 6);
    end

    m1=  [ image_max_peak_mean{1}{1}(curr_roi,3)  ...
        image_max_peak_mean{1}{2}(curr_roi,3) ...
        image_max_peak_mean{2}{1}(curr_roi,2)  ...
        image_max_peak_mean{2}{2}(curr_roi,2)];

    e1=[ image_max_peak_error{1}{1}(curr_roi,3)  ...
        image_max_peak_error{1}{2}(curr_roi,3)...
        image_max_peak_error{2}{1}(curr_roi,2)  ...
        image_max_peak_error{2}{2}(curr_roi,2)];

    image_each_mice=  { permute(image_max_peak{1}{1}(curr_roi,3,:),[3,2,1])  ...
        permute(image_max_peak{1}{2}(curr_roi,3,:),[3,2,1])  ...
        permute( image_max_peak{2}{1}(curr_roi,2,:),[3,2,1])   ...
        permute(image_max_peak{2}{2}(curr_roi,2,:),[3,2,1]) };

    p_val{curr_roi}=arrayfun(@(id) ds.shuffle_test( image_each_mice{id},image_each_mice{id+1},0,2),[1 3],'UniformOutput',true)

    hold on
    for curr_draw=1:2
        errorbar(a2, [curr_draw curr_draw+2], m1([curr_draw curr_draw+2]), e1([curr_draw curr_draw+2]), ...
            'o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
            'MarkerFaceColor', colors{curr_draw}, 'Color', colors{curr_draw}, ...
            'LineStyle','none','CapSize',0);
    end
    title(roi1(curr_roi).name,'FontWeight','normal')

    set(a2,'XLim',[0.5 4.5],'YLim',[0 scale_plot], ...
        'XTick',[1.5 3.5],'XTickLabel',{'V passive','A passive'}, ...
        'Box','off','Color','none');
    ylabel(a2,'\Delta F/F_{0}')

    % inset：左上角显示 ROI mask
    mp = a2.Position;

    insetAx = axes('Position', [mp(1)+0.050*mp(3), mp(2)+0.75*mp(4), 0.30*mp(3), 0.30*mp(4)]);
    imagesc(insetAx, roi1(curr_roi).data.mask); ap.wf_draw('ccf',[0.5 0.5 0.5]);
    axis(insetAx,'image','off'); ylim(insetAx,[0 200]); xlim(insetAx,[20 220]);
    clim(insetAx,[0 1]); colormap(insetAx, ap.colormap('WK')); uistack(insetAx,'bottom');
end


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s4d.eps'), ...
    'ContentType','vector');
 clearvars('-except',main_preload_vars{:});
