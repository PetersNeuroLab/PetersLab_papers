%% fig 3 context difference
main_preload_vars = who;
load(fullfile(Path,'data','wf_passive_kernels'));
load(fullfile(Path,'data','wf_task_kernels'));
load(fullfile(Path,'data','behavior'));


temp_plot_task=cellfun(@(x) permute(ds.make_each_roi( cat(4,nanmean(plab.wf.svd2px(U_master(:,:,1:170),x(:,:,1:3,:) ),4),...
    nanmean(plab.wf.svd2px(U_master(:,:,1:170),x(:,:,4:8,:) ),4)),...
    length(t_kernels),roi1([1 3])),[1,2,4,3]),...
    wf_task_kernels_across_day,'UniformOutput',false)
temp_plot_passive=cellfun(@(x,id,stim) ...
    permute(ds.make_each_roi( cat(5,nanmean(plab.wf.svd2px(U_master(:,:,1:size(x{id},1)),x{id}(:,:,stim,4:6,:) ),5),...
    nanmean(plab.wf.svd2px(U_master(:,:,1:size(x{id},1)),x{id}(:,:,stim,7:11,:) ),5)),...
    length(t_kernels),roi1([1 3])),[1,2,5,4,3]),...
    wf_passive_kernels_across_day,{1;2},{3;2},'UniformOutput',false);



temp_task=cellfun(@(x) ds.make_each_roi( nanmean(plab.wf.svd2px(U_master(:,:,1:170),x(:,:,1:8,:) ),5),length(t_kernels),roi1),...
    wf_task_kernels_across_day,'UniformOutput',false);
temp_passive=cellfun(@(x,id,stim) ...
    permute(ds.make_each_roi( nanmean(plab.wf.svd2px(U_master(:,:,1:size(x{id},1)),x{id}(:,:,stim,4:11,:) ),6),length(t_kernels),roi1),[1,2,4,3]),...
    wf_passive_kernels_across_day,{1;2},{3;2},'UniformOutput',false);

workflow={'visual position';'audio volume'};

temp_task_roi =cellfun(@(k,w,name) cellfun(@(k1,w1)  cellfun(@(x) ...
    ds.make_each_roi( plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),length(t_kernels),roi1),k1(ismember(w1,name)),'UniformOutput',false)...
    ,k,w.workflow_name,'UniformOutput',false),wf_task_kernel_each_mice,behavior_each_mice,workflow,'UniformOutput',false);
temp_task_peak=cellfun(@(a1) cellfun(@(a2)  cellfun(@(a3)  max(a3([1 3],kernels_period),[],2),a2,...
    'UniformOutput',false),a1, 'UniformOutput',false),temp_task_roi,'UniformOutput',false);
temp_task_peak2=cellfun(@(a1) cellfun(@(a2)  cat(2,a2{:})'  ,a1, 'UniformOutput',false),temp_task_peak,'UniformOutput',false);


temp_passive_roi =cellfun(@(k,w,name,id) cellfun(@(k1,w1)  cellfun(@(x) ...
    ds.make_each_roi( plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),length(t_kernels),roi1),k1(find(ismember(w1,name))+3),'UniformOutput',false)...
    ,k{id},w.workflow_name,'UniformOutput',false),wf_passive_kernel_each_mice,behavior_each_mice,workflow,{1;2},'UniformOutput',false);


temp_passive_peak=cellfun(@(a1,stim) cellfun(@(a2)  cellfun(@(a3)  max(a3([1 3],kernels_period,stim),[],2),a2,...
    'UniformOutput',false),a1, 'UniformOutput',false),temp_passive_roi,{3;2},'UniformOutput',false);


temp_passive_peak2=cellfun(@(a1) cellfun(@(a2)  cat(2,a2{:})'  ,a1, 'UniformOutput',false),temp_passive_peak,'UniformOutput',false);
temp_passive_peak3=cellfun(@(a1)   cat(1,a1{:}) ,temp_passive_peak2,'UniformOutput',false);

temp_perform =cellfun(@(w,name) cellfun(@(p1,w1) p1(ismember(w1,name),1),...
    w.performance ,w.workflow_name,'UniformOutput',false),behavior_each_mice,workflow,'UniformOutput',false);

temp_learn =cellfun(@(w,name) cellfun(@(p1,w1) p1(ismember(w1,name),1),...
    w.learned ,w.workflow_name,'UniformOutput',false),behavior_each_mice,workflow,'UniformOutput',false);

%
used_area=[1  3 ];
title_images={'pre','post'};
title_area={'l-mPFC','r-mPFC','l-aPFC','r-aPFC','l-PPC','r-PPC','all-PFC','auditory area','','','','V1'}
scale1=0.0004;
%
Color={'B','R'};
colors{1} = { [0 0 1],[1 0 0]}; % 定义颜色
colors{2} = { [0.5 0.5 1],[1 0.5 0.5]}; % 定义颜色
color_roi={[1 0.5 1],[ 0.5 1 0.5]}
figure('Position',[50 50 1800 500])
t1 = tiledlayout(length(used_area), 8, 'TileSpacing', 'tight', 'Padding', 'tight');
font_size=12

for curr_roi=1:2
    curr_area = used_area(curr_roi);  % 显式编号
    a1=nexttile(t1,curr_roi*8-7)
    imagesc(roi1(curr_area).data.mask )
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    axis image off

    ylim([0 200])
    xlim([20 220])
    clim( [ 0, 1]);
    colormap( a1,ap.colormap('WK'));
    title(roi1(curr_area).name,'FontWeight','normal', 'Interpreter', 'none')
    a1.FontSize = font_size;

end

for curr_roi=1:2
    for curr_group=1:2

        use_roi=used_area(curr_roi);
        a1=nexttile(t1,curr_roi*8-8+3*curr_group-1)
        imagesc(t_kernels,[],permute(temp_task{curr_group}(use_roi,:,:),[2,3,1])');
        colormap(a1, ap.colormap(['W' Color{curr_group}]));
        ylim([0.5 8.5])
        xlim([-0.2 0.5])
        yline(3.5);
        clim(scale1 .* [0, 1]);
        xticks([-0.2 0 0.2 0.5])

        if curr_roi==1

        else
            xlabel('time (s)')
        end


        yticks([2  6 ]); % 设置 y 轴的刻度位置（2代表naive stage中间位置，8代表stage1中间位置）
        yticklabels(title_images); % 设置对应的标签
        % ytickangle(90);
        ylabel('days')% 旋转 90 度

        if curr_roi==1
            title('task','FontWeight','normal','FontSize',20)
        end
        a1.FontSize = font_size;

        a1=nexttile(t1,curr_roi*8-8+3*curr_group)
        imagesc(t_kernels,[],permute(temp_passive{curr_group}(use_roi,:,:),[2,3,1])');
        colormap(a1, ap.colormap(['W' Color{curr_group}]));
        ylim([0.5 8.5])
        xlim([-0.2 0.5])
        yline(3.5);
        clim(scale1 .* [0, 1]);
        yticks([])
        xticks([-0.2 0 0.2 0.5])

        if curr_roi==1

            cb = colorbar('southoutside');  % 横向放在下方
            pos = cb.Position;   % [left bottom width height]
            pos(4) = pos(4)/2;   % 缩短高度
            pos(3) = pos(3) /2;   % 缩短高度
            pos(1) =pos(1)+0.05;   % 缩短高度
            pos(2) =pos(2)-0.01;   % 缩短高度

            cb.Position = pos;
            cb.Ticks = [cb.Limits(1), cb.Limits(2)];   % 只显示最小和最大
            cb.Label.String = '\DeltaFR/FR_{0}';   % 给 colorbar 加标签

        else
            xlabel('time (s)')

        end
        if curr_roi==1
            title('passive','FontWeight','normal','FontSize',14)
        end
        a1.FontSize = font_size;

    end
end

colors_1{1}=[[0 0 1];[0.5 0.5 1]];
colors_1{2}=[[1 0 0];[1 0.5 0.5]];
slope_task=cell(2,1);
slope_pass=cell(2,1);
a1=cell(2,2)
for curr_roi=1:2
    for curr_group=1:2
        a1{curr_roi,curr_group}=nexttile(t1,curr_roi*8-8+curr_group*3+1)
        hold on
        h1 = scatter(NaN,NaN,20,'filled','MarkerFaceColor',[0.2 0.2 0.2],'LineWidth',1);

        cellfun(@(p2,t2,l2)   scatter(p2(l2==0),t2(l2==0,curr_roi),20,'filled',...
            'MarkerFaceColor',[0.2 0.2 0.2],'LineWidth',1),...
            temp_perform{curr_group},temp_task_peak2{curr_group},temp_learn{curr_group},'UniformOutput',false )
        h2 = scatter(NaN,NaN,20,'filled','MarkerFaceColor',colors_1{curr_group}(1,:),'LineWidth',1);
        cellfun(@(p2,t2,l2)   scatter(p2(l2==1),t2(l2==1,curr_roi),20,'filled',...
            'MarkerFaceColor',colors_1{curr_group}(1,:),'LineWidth',1),...
            temp_perform{curr_group},temp_task_peak2{curr_group},temp_learn{curr_group},'UniformOutput',false )


        h3 = scatter(NaN,NaN,20,'filled','MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1);
        cellfun(@(p2,t2,l2)   scatter(p2(l2==0),t2(l2==0,curr_roi),20,'filled',...
            'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1),...
            temp_perform{curr_group},temp_passive_peak2{curr_group},temp_learn{curr_group},'UniformOutput',false )
        h4 = scatter(NaN,NaN,20,'filled','MarkerFaceColor',colors_1{curr_group}(2,:),'LineWidth',1);
        cellfun(@(p2,t2,l2)   scatter(p2(l2==1),t2(l2==1,curr_roi),20,'filled',...
            'MarkerFaceColor',colors_1{curr_group}(2,:),'LineWidth',1),...
            temp_perform{curr_group},temp_passive_peak2{curr_group},temp_learn{curr_group},'UniformOutput',false )



        learn_3=logical(cat(1,temp_learn{curr_group}{:}));

        task_peak3=cat(1,temp_task_peak2{curr_group}{:});
        perform3=cat(1,temp_perform{curr_group}{:});

        p_task = polyfit(perform3(learn_3), task_peak3(learn_3,curr_roi), 1);
        x_fit_task = linspace(0, 1, 2);
        y_fit_task = polyval(p_task, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors_1{curr_group}(1,:));


        passive_peak3=cat(1,temp_passive_peak2{curr_group}{:});


        p_passive = polyfit(perform3(learn_3), passive_peak3(learn_3,curr_roi), 1);
        x_fit_passive = linspace(0, 1, 2);
        y_fit_passive = polyval(p_passive, x_fit_passive);
        plot(x_fit_passive, y_fit_passive, '-', 'LineWidth', 2,'Color',colors_1{curr_group}(2,:));


        [R_task,P_task] = corr(perform3(learn_3), task_peak3(learn_3,curr_roi));

        [R_passive,P_passive] = corr(perform3(learn_3),  passive_peak3(learn_3,curr_roi));


        slope_task{curr_group}{curr_roi}= cellfun(@(perform,peak,learned) diff(polyval( polyfit( perform, peak(:,curr_roi),1), linspace(0, 1, 2))),...
            temp_perform{curr_group},temp_task_peak2{curr_group},temp_learn{curr_group},'UniformOutput',true);

        slope_pass{curr_group}{curr_roi}= cellfun(@(perform,peak,learned) diff(polyval( polyfit( perform, peak(:,curr_roi),1), linspace(0, 1, 2))),...
            temp_perform{curr_group},temp_passive_peak2{curr_group},temp_learn{curr_group},'UniformOutput',true);


        ylim([0 0.0005])
        xlim([-0.1 1])
        % title(roi_name{curr_roi} ,'FontWeight','normal')
        ylabel('max \Delta F/F_{0}')

        xlabel('performance')
        axis square
        a1{curr_roi,curr_group}.FontSize = font_size;
        if curr_roi==1
            legend([h1 h2 h3 h4], ...
                {'task pre','task post','passive pre','passive post'},'NumColumns',2, ...
                'Location','northoutside','Box','off');
        end
        set(gca, 'Color', 'none');        % 坐标轴背景透明

    end
end

for curr_roi=1:2
    for curr_group=1:2

        mainPos = get(a1{curr_roi,curr_group}, 'Position');  % [left bottom width height]

        % 计算 inset 的位置（嵌在当前 tile 的左上角）
        inset_width = 0.4* mainPos(3);    % inset 占 tile 宽度的 30%
        inset_height = 0.4 * mainPos(4);   % inset 占 tile 高度的 30%
        inset_left = mainPos(1) - 0* mainPos(3);  % tile 左侧偏右一点
        inset_bottom = mainPos(2) + 0.65 * mainPos(4); % tile 底部偏上
        insetAx = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
        ap.errorfill(t_kernels, nanmean(temp_plot_task{curr_group}(curr_roi,:,:,1),3),...
            std(temp_plot_task{curr_group}(curr_roi,:,:),0,3)./sqrt(size(temp_plot_task{curr_group},3)),...
            [0 0 0],0.1,0.1 )
        ap.errorfill(t_kernels, nanmean(temp_plot_passive{curr_group}(curr_roi,:,:,1),3),...
            std(temp_plot_passive{curr_group}(curr_roi,:,:),0,3)./sqrt(size(temp_plot_passive{curr_group},3)),...
            [0.5 0.5 0.5],0.1,1 )

        ap.errorfill(t_kernels, nanmean(temp_plot_task{curr_group}(curr_roi,:,:,2),3),...
            std(temp_plot_task{curr_group}(curr_roi,:,:),0,3)./sqrt(size(temp_plot_task{curr_group},3)),colors_1{curr_group}(1,:),0.1,0.1 )
        ap.errorfill(t_kernels, nanmean(temp_plot_passive{curr_group}(curr_roi,:,:,2),3),...
            std(temp_plot_passive{curr_group}(curr_roi,:,:),0,3)./sqrt(size(temp_plot_passive{curr_group},3)),colors_1{curr_group}(2,:),0.1,1 )
        xlim([-0.1 0.5])
        ylim([-1e-4 3e-4])
        axis off

        uistack(insetAx, 'bottom');
    end
end



colors1={ [0 0 1];[0.5 0.5 1];[1 0 0];[1 0.5 0.5]};
for curr_roi=1:2
    a2=nexttile(t1,curr_roi*8)


    temp_plot_data=[cellfun(@(x)  x(curr_roi)  , slope_task,'UniformOutput',true );...
        cellfun(@(x)  x(curr_roi)  , slope_pass,'UniformOutput',true )];

    ds.make_bar_plot(temp_plot_data([1 3 2 4]),'ColorCell',colors1,'BarAlpha',1,'ShowDots',0,'CentralTendency','median','ShowErrorCaps',0)
    ds.make_bar_plot(temp_plot_data([1 3 2 4]),'ColorCell',colors1,'BarAlpha',1,'ShowDots',0,'CentralTendency','median','ShowErrorCaps',0)

    hold on

    xline(2.5,'LineStyle',':')
    xticks([1:4])
    xticklabels({'V-task','V-passive','A-task','A-passive'})

    if curr_roi==1
        y_sig = 0.00025;

        ylim([-0.0001 0.0003])

    else
        y_sig = 0.0006;

        ylim([-0.0001 0.00065])
    end

    ylabel('slope')
    set(gca, 'Color', 'none');        % 坐标轴背景透明

    p_in_group=cellfun(@(a,b) ds.shuffle_test(a{curr_roi},b{curr_roi},0,1),slope_task,slope_pass);
    p_across_group_task= ds.shuffle_test(slope_task{1}{curr_roi},slope_task{2}{curr_roi},0,1);
    p_across_group_pass= ds.shuffle_test(slope_pass{1}{curr_roi},slope_pass{2}{curr_roi},0,1);
    temp_p=[p_in_group;p_across_group_task;p_across_group_pass];
    lines={[1 2],[3 4],[1 3],[2 4]};

    temp_gap=0.02*y_sig;
    for curr_i=1:4
        if temp_p(curr_i) < 0.05
            % stars = repmat('*',1,sum(temp_p(curr_i)<[0.05 0.01 0.001]));
            stars = '*';

            plot(lines{curr_i}, [1 1]*(y_sig+temp_gap*curr_i), 'k-');
            text(mean(lines{curr_i}), (y_sig+temp_gap*curr_i)+temp_gap, stars, 'HorizontalAlignment','center');
        end
        drawnow

    end


    a2.FontSize = 12;

end


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 3.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
clearvars('-except',main_preload_vars{:});
