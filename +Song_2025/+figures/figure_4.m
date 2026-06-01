%% fig 4  task kernels cross modality

main_preload_vars = who;

load_dataset='wf_task_kernels';
load(fullfile(Path,'data',load_dataset));
load(fullfile(Path,'data','behavior'));
% behavior

reaction_time_1=cellfun(@(x)  x.reaction_time([1:8 15:20],:) ,behavior_across_day,'UniformOutput',false)
reaction_time_mean=cellfun(@(x) nanmean( x.reaction_time([1:8 15:20],:),2) ,behavior_across_day,'UniformOutput',false)
reaction_time_error=cellfun(@(x) std( x.reaction_time([1:8 15:20],:),0,2,'omitmissing')./sqrt(size(x.reaction_time([1:8 15:20],:),2)),...
    behavior_across_day,'UniformOutput',false);


% kernels
tem_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[ 7 8 13 14],:)),  wf_task_kernels_across_day,'UniformOutput',false);

image_max_eachmice=cellfun(@(x)    cat(3,  nanmean(max(x(:,:,kernels_period,[1 2],:),[],3),[4 5]),...
    nanmean(max(x(:,:,kernels_period,[3 4],:),[],3),[4 5])),tem_image,'UniformOutput',false  );


tem_image_across_day=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[1:8 15:20],:)),  wf_task_kernels_across_day,'UniformOutput',false);

tem_roi_across_day= cellfun(@(x) ds.make_each_roi(permute(max(x(:,:,kernels_period,:,:),[],3),[1,2,4,5,3]) ,size(x,3),roi1),...
    tem_image_across_day,'UniformOutput',false)
tem_roi_across_day_mean=cellfun(@(x) nanmean(x,3) ,tem_roi_across_day,'UniformOutput',false );
tem_roi_across_day_error=cellfun(@(x) std(x,0,3,'omitmissing')./sqrt(size(x,3)) ,tem_roi_across_day,'UniformOutput',false );


%

scale_image=0.0004;
scale_plot=0.0005;
Color={'G','P'};
colors = [ ...
    84 130 53  % #548235
    112  48 160  % #7030A0
    ] / 255;

figure('Position', [50 50 450 450]);
t = tiledlayout(3, 3, 'TileSpacing', 'tight', 'Padding', 'none');
order={[1 2],[2 1]}
a1=cell(2,2)
for curr_image=1:2

    for curr_group=1:2

        a1{curr_image,curr_group}=nexttile(curr_group+3*curr_image)
        imagesc(image_max_eachmice{curr_group}(:,:,order{curr_group}(curr_image)))
        axis image off;
        clim(scale_image .* [0.25, 0.75]);
        colormap(a1{curr_image,curr_group}, ap.colormap(['W' Color{curr_group}] ));
        clim(scale_image .* [0, 1]);
        ap.wf_draw('ccf', [0.5 0.5 0.5]);
        % title(titles{curr_image},'FontWeight','normal')

    end
end

for curr_group=1:2
    cb = colorbar(a1{1,curr_group},'southoutside');  % 横向放在下方
    pos = cb.Position;   % [left bottom width height]
    pos(4) = pos(4)/2;   % 缩短高度
    pos(3) = pos(3) /2;   % 缩短高度
    % pos(2) =0.04;   % 缩短高度

    cb.Position = pos;
    cb.Ticks = [cb.Limits(1), cb.Limits(2)];   % 只显示最小和最大

    cb.Label.String = '\DeltaFR/FR_{0}';   % 给 colorbar 加标签

end



kernels_p=cell(3,1);
% use_area=[1 3];
for curr_area=[1 3]
    ax1=nexttile(t,1.5*curr_area+4.5)
    style={'-','--'}
    for curr_group=1:2
        ap.errorfill(1:8,tem_roi_across_day_mean{curr_group}(curr_area,1:8),tem_roi_across_day_error{curr_group}(curr_area,1:8) ,colors(curr_group,:),0.1,0);
        ap.errorfill(9:14,tem_roi_across_day_mean{curr_group}(curr_area,9:14),tem_roi_across_day_error{curr_group}(curr_area,9:14) ,colors(curr_group,:),0.1,0);
        plot(1:8,tem_roi_across_day_mean{curr_group}(curr_area,1:8),'Color', colors(curr_group,:),'LineStyle',style{curr_group},'LineWidth',2);
        plot(9:14,tem_roi_across_day_mean{curr_group}(curr_area,9:14),'Color',colors(curr_group,:),'LineStyle',style{3-curr_group},'LineWidth',2);

        xlim(ax1,[3 14])
        xticks(ax1,[3 8 9 14 ]); % 设置 y 轴的刻度位置（2代表naive stage中间位置，8代表stage1中间位置）
        xticklabels(ax1,{'-6','-1','0','4'}); % 设置对应的标签
        ylim(ax1,[0 scale_plot])
        ylabel(ax1,'\Delta F/F_{0}')
        xlabel('day from transfer')

    end
    title(roi1(curr_area).name,'FontWeight','normal')
    kernels_p{curr_area}=cellfun(@(x) ds.shuffle_test(permute(nanmean(x(curr_area,7:8,:),2),[1,3,2]),...
        permute(nanmean(x(curr_area,9:10,:),2),[1,3,2]),1,2),tem_roi_across_day,'UniformOutput',true);

    test_p=cellfun(@(x) ds.shuffle_test(permute(nanmean(x(curr_area,7:8,:),2),[1,3,2]),...
        permute(nanmean(x(curr_area,9:10,:),2),[1,3,2]),1,2),tem_roi_across_day,'UniformOutput',true)>0.95


    % test_p 是 [val1 val2]
    offset = 0;  % 用于竖直堆叠
    y_base=[0.0002 0 0.0004]
    if test_p(1)
        text(9.5, y_base(curr_area) + offset, '*', 'Color',colors(1,:), 'FontSize',14, ...
            'HorizontalAlignment','center');
        offset = offset + 0.0001;  % 往上移一格
    end
    if test_p(2)
        text(9.5, y_base(curr_area)+ offset, '*', 'Color',colors(2,:), 'FontSize',14, ...
            'HorizontalAlignment','center');
    end
    xline(8.5,'LineStyle',':','LineWidth',1,'Color',[0.5 0.5 0.5])



    set(gca,'Color','none')

    mainPos = get(ax1, 'Position');  % [left bottom width height]
    % 计算 inset 的位置（嵌在当前 tile 的左上角）
    inset_width = 0.3 * mainPos(3);    % inset 占 tile 宽度的 30%
    inset_height = 0.3 * mainPos(4);   % inset 占 tile 高度的 30%
    inset_left = mainPos(1) + 0.05 * mainPos(3);  % tile 左侧偏右一点
    inset_bottom = mainPos(2) + 0.65 * mainPos(4); % tile 底部偏上
    insetAx = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
    imagesc(roi1(curr_area).data.mask )
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    axis image off
    ylim([0 200])
    xlim([20 220])
    clim( [ 0, 1]);
    colormap( insetAx,ap.colormap('WK'));
    uistack(insetAx, 'bottom');

end


a1=nexttile(t,3)
for curr_group=1:2
    hold on
    ap.errorfill(1:8, reaction_time_mean{curr_group}(1:8),...
        reaction_time_error{curr_group}(1:8),colors(curr_group,:),0.1,0)
    ap.errorfill(9:14, reaction_time_mean{curr_group}(9:14),...
        reaction_time_error{curr_group}(9:14),colors(curr_group,:),0.1,0)

    h(1)=plot(1:8,reaction_time_mean{curr_group}(1:8),'Color', colors(curr_group,:),'LineStyle',style{curr_group},'LineWidth',2);
    h(2)= plot(9:14,reaction_time_mean{curr_group}(9:14),'Color',colors(curr_group,:),'LineStyle',style{3-curr_group},'LineWidth',2);


    set(gca,'Color','none')

end
set(gca, 'YScale', 'log');   % 把 y 轴改成对数刻度
xline(8.5,'LineStyle',':','LineWidth',1,'Color',[0.5 0.5 0.5])
xlim([3 14])
xticks([3 8 9 14 ]); % 设置 y 轴的刻度位置（2代表naive stage中间位置，8代表stage1中间位置）
xticklabels({'-6','-1','0','4'}); % 设置对应的标签
ylabel('Reaction time (s)')
legend([h(1), h(2)], {'auditory', 'visual'} ,'Location','northwest','Box','off');

% legend({'',legned_name{1},'','','','','',legned_name{2}},'Location','northoutside','Box','off','Orientation','horizontal')
ylim([0 1])
yticks([0.1 0.2 0.4 0.6 0.8 1])
% yticklabels({'0.1','0.32','1'})
xlabel('day from transfer')
p_val_corss_group=cellfun(@(x)  ds.shuffle_test (nanmean(x(7:8,:),1),nanmean(x(9:10,:),1),1,2 )  , reaction_time_1,'UniformOutput',true)
test_p= cellfun(@(x)  ds.shuffle_test (nanmean(x(7:8,:),1),nanmean(x(9:10,:),1),1,2 )>0.95   , reaction_time_1,'UniformOutput',true)
% test_p 是 [val1 val2]
offset = 0;  % 用于竖直堆叠
y_base = 0.75;   % 星号基准高度，可以根据数据调节

if test_p(1)
    text(9.5, y_base + offset, '*', 'Color',colors(1,:), 'FontSize',14, ...
        'HorizontalAlignment','center');
    offset = offset + y_step;  % 往上移一格
end
if test_p(2)
    text(9.5, y_base + offset, '*', 'Color',colors(2,:), 'FontSize',14, ...
        'HorizontalAlignment','center');
end
set(gca,'Color','none')
exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 4.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
clearvars('-except',main_preload_vars{:});
