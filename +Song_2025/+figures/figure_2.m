%% fig 2 a-d  passive kernels images
main_preload_vars = who;
load_dataset='wf_passive_kernels';
load(fullfile(Path,'data',load_dataset));


oder={[3 1 2],[2 1 3]};

scale=0.0003;
xlabel_all={'left','center','right';'4k Hz','8k Hz','12k Hz'}
line_colors={[0.5 0.5 0.5],[0.5 0.5 0.5],[0 0 1];[0.5 0.5 0.5],[1 0 0],[0.5 0.5 0.5]}
figure('Position',[50 50 450 300])
mainLayout = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'none');

for curr_group=1:2
    for  workflow_idx =curr_group
        main_preload_vars = who;
        imageLayout = tiledlayout(mainLayout, 1, 4, ...
            'TileSpacing', 'tight', 'Padding', 'tight');
        imageLayout.Layout.Tile = curr_group;  % 明确放在主 layout 的第 1 个 tile
        used_area=workflow_idx*2-1;
        used_stim=4-workflow_idx;
        scale=0.0003;
        Color={'B','R'};
        % darw iamge
        temp_image= plab.wf.svd2px(U_master(:,:,1:size(wf_passive_kernels_across_day{workflow_idx}{curr_group},1)),wf_passive_kernels_across_day{workflow_idx}{curr_group}(:,:,:,10:11,:));
        buf_images_1= permute(nanmean(max(temp_image(:,:,kernels_period,:,:,:),[],3),[5 6]),[1,2,4,3,5 6]);

        for curr_stim=1:3
            use_stim=oder{curr_group}(curr_stim);
            ax(curr_group,curr_stim)=nexttile(imageLayout);

            imagesc(buf_images_1(:,:,use_stim))
            axis image off;
            ap.wf_draw('ccf', [0.5 0.5 0.5]);
            clim(scale .* [ 0, 1]);
            colormap(ax(curr_group,curr_stim), ap.colormap(['W' Color{workflow_idx}]));
            title(xlabel_all{workflow_idx,use_stim},'FontWeight','normal','FontSize',10)

        end

    end
end

for curr_group=1:2
    cb=colorbar(ax(curr_group,curr_stim));
    cb.Units = 'normalized';
    cb.Position = [0.8, 1-curr_group*0.4 ,0.02, 0.1];
    % cb.Position
end


% ap.prettyfig

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 2a_d.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
clearvars('-except',main_preload_vars{:});


%% fig 2 ef  passive kernels selectivity
main_preload_vars = who;
load_dataset='wf_passive_kernels';
load(fullfile(Path,'data',load_dataset));
line_colors={[0.7 0.7 1],[0.7 0.7 1],[0 0 1];[1 0.7 0.7],[1 0 0],[1 0.7 0.7]};


temp_roi_plot_mean=cell(2,1);
temp_roi_plot_error=cell(2,1);
temp_roi_peak_mean=cell(2,1);
temp_roi_peak_error=cell(2,1);
temp_roi_peak=cell(2,1);
for curr_group=1:2
    for  workflow_idx =curr_group

        temp_image= plab.wf.svd2px(U_master(:,:,1:size(wf_passive_kernels_across_day{workflow_idx}{curr_group},1)),...
            wf_passive_kernels_across_day{workflow_idx}{curr_group}(:,:,:,10:11,:));

        temp_each_roi=ds.make_each_roi(temp_image, length(t_kernels),roi1);
        temp_roi_plot_mean{curr_group}=nanmean(temp_each_roi,[4,5]);
        temp_roi_plot_error{curr_group}=std(nanmean(temp_each_roi,4),0,5)./sqrt(size(temp_each_roi,5));
        temp_roi_peak{curr_group}=permute(nanmean(max(temp_each_roi(:,kernels_period,:,:,:),[],2),4),[1,3,5,2,4]);

        temp_roi_peak_mean{curr_group}=permute(nanmean(max(temp_each_roi(:,kernels_period,:,:,:),[],2),[4,5]),[1,3,2]);
        temp_roi_peak_error{curr_group}=permute(std(nanmean(max(temp_each_roi(:,kernels_period,:,:,:),[],2),4),0,5)./sqrt(size(temp_each_roi,5)),[1,3,2]);

    end
end

figure('Position',[50 50 900 400])
plot_layout=tiledlayout(1,4, 'TileIndexing','columnmajor','TileSpacing','tight','Padding','tight')

use_area=[1 3];
colors={[0.7 0.7 1],[1 0.7 0.7]}
oder={[3 1 2],[2 1 3]};

for area_idx = 1:2  % 对应 curr_area = [1 3]
    sub_fig=tiledlayout(plot_layout,3,1,'TileSpacing','none','Padding','none');
    sub_fig.Layout.Tile=2*area_idx-1;

    curr_area = use_area(area_idx);  % 显式编号
    col_idx = area_idx;  % 当前是在第几列放置（因为 columnmajor）
    tt=nexttile(sub_fig);
    imagesc(roi1(curr_area).data.mask );
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    axis image off

    ylim([0 200]);
    xlim([20 220]);
    clim( [ 0, 1]);
    colormap( tt,ap.colormap('WK'));
    title(roi1(curr_area).name,'FontWeight','normal');

    for curr_group=1:2

        used_stim=4-curr_group;
        nexttile(sub_fig);
        hold on
        for curr_stim=1:3
            ap.errorfill(t_kernels,temp_roi_plot_mean{curr_group}(curr_area,:,curr_stim),...
                temp_roi_plot_error{curr_group}(curr_area,:,curr_stim),line_colors{curr_group,curr_stim},0.1,1,2);
            xlim([-0.05 0.4]);
            ylim(1e-4*[-0.3 2]);
        end

        axis off
    end

    line_loc=-0.3*1e-4;
    line([-0.05 0],[line_loc line_loc],'Color',[0 0 0],'LineStyle','-');
    line([-0.05 -0.05],[line_loc line_loc+0.5e-4],'Color',[0 0 0],'LineStyle','-');
    text(0.15, -0.1,'0.05s' ,  'Units','normalized', 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'FontSize',10, 'FontWeight','normal');
    text(-0.25, 0.25,{'0.5\times10^{-4}', '\DeltaF/F_{0}'} ,  'Units','normalized', 'HorizontalAlignment','center', ...
        'VerticalAlignment','top', 'FontSize',10, 'FontWeight','normal', 'Rotation',90);



    nexttile(plot_layout,area_idx*2);  % 在第2块区域的后面画一整栏（例如列3、列4）

    % oder{curr_group}(curr_stim)

    hold on
    for curr_group=1:2
        plot(1:3,temp_roi_peak_mean{curr_group}(curr_area,oder{curr_group}),'Color',colors{curr_group}, 'LineWidth', 2)

        for  curr_stim=1:3
            errorbar( curr_stim ,temp_roi_peak_mean{curr_group}(curr_area,oder{curr_group}(curr_stim)),...
                temp_roi_peak_error{curr_group}(curr_area,oder{curr_group}(curr_stim))...
                ,'-o','MarkerSize',4, 'LineWidth', 2,'Color',line_colors{curr_group,oder{curr_group}(curr_stim)},'MarkerFaceColor',line_colors{curr_group,oder{curr_group}(curr_stim)},'CapSize',0)
        end
    end
    xticks([1:3]); % 设置 y 轴的刻度位置（2代表naive stage中间位置，8代表stage1中间位置）
    xticklabels({'R/8K','L/4K','C/12K'}); % 设置对应的标签
    xlabel('d')
    xlim([0.5 3.5])
    scale=0.00025;

    ylim(scale .* [-0, 1 ]);
    % title(roi1(curr_area).name,'FontWeight','normal','FontSize',10)
    yticks(1e-4*[0 1 2])

    ylabel('\DeltaF/F_{0}')
    % xlabel('stim types')
    box off
    set(gca,'Color','none')


end


use_area=[1 3];
p_vals=nan(2,2);
for curr_group=1:2
    switch curr_group
        case 1
            idx={3 ,[1 2]}
        case 2
            idx={2, [1 3]}
    end
    for area_idx = 1:2
        curr_area=use_area(area_idx)
        temp_dat=permute(temp_roi_peak{curr_group}(curr_area,:,:),[3,2,1]);
        p_vals(curr_group,area_idx) = ds.shuffle_test(temp_dat(:,idx{1}) ,nanmean(temp_dat(:,idx{2}),2),0,1)
    end
end

ap.prettyfig

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 2e_f.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
clearvars('-except',main_preload_vars{:});