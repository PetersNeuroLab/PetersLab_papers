%% Generate figures for Song et al 2025
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



%%  fig 1c_d  task kernels images
main_preload_vars = who;

load_dataset='wf_task_kernels';
load(fullfile(Path,'data',load_dataset));
tem_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[2 3 7 8],:)),  wf_task_kernels_across_day,'UniformOutput',false);
image_max_eachmice=cellfun(@(x)    cat(3,  nanmean(max(x(:,:,kernels_period,[1 2],:),[],3),[4 5]),...
    nanmean(max(x(:,:,kernels_period,[3 4],:),[],3),[4 5])),tem_image,'UniformOutput',false  );


% drawnow
figure('Position',[50 50 700 600])
imagelayout = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

colors = {'B','R'};
titles = {'pre learned','well trained'};

axs = gobjects(2,2);                 % 存轴句柄
for curr_group=1:2
    for curr_stage=1:2
        ax = nexttile(imagelayout);
        axs(curr_group,curr_stage) = ax;
        imagesc(ax, image_max_eachmice{curr_group}(:,:,curr_stage))
        axis(ax,'image','off')
        clim(ax, 0.0003*[0,1]);
        ap.wf_draw('ccf',[0.5 0.5 0.5]);
        colormap(ax, ap.colormap(['W' colors{curr_group}]));
        if curr_group==1
            title(ax, titles{curr_stage},'FontSize',10,'FontWeight','normal')
        end
    end
end

% —— 每一行右侧放一个 colorbar（图像外面），高度=子图高度的1/3 ——
drawnow;                             % 先让布局稳定
gap = 0.01;                          % 与子图的水平间距
cbw = 0.02;                          % colorbar 的宽度（归一化坐标）
shrink = 1/3;                        % 高度比例

for r = 1:2
    ax = axs(r,2);                   % 每行最后一个（右侧）子图
    p = ax.Position;                 % [x y w h] 归一化
    cb = colorbar(ax);               % 关联该轴（继承其 CLim 和 colormap）
    cb.Units = 'normalized';
    h = p(4)*shrink;                 % colorbar 高度=子图高度的1/3
    x = p(1) + p(3) + gap-0.01;           % 放在子图右侧，留一点间隙
    y = p(2) ;           % 垂直居中（也可改为 y=p(2) 贴底）
    cb.Position = [x, y, cbw, h];
end
drawnow;

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 1 cd.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});

%% fig 1e  task difference with permutation test
main_preload_vars = who;

load_dataset='wf_task_kernels';
load(fullfile(Path,'data',load_dataset));

tem_image_2=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[7 8],:)),  wf_task_kernels_across_day,'UniformOutput',false);
image_max_eachmice=cellfun(@(x)   permute(nanmean(max(x(:,:,kernels_period,[1 2],:),[],3),4),[1,2,5,3,4]),...
    tem_image_2,'UniformOutput',false  );

image_max_base_eachmice=cellfun(@(x)   permute(nanmean(max(x(:,:,t_kernels<0 &t_kernels>-0.3,[1 2],:),[],3),4),[1,2,5,3,4]),...
    tem_image_2,'UniformOutput',false  );


% permutation test
threshold=0.0001;
A=image_max_base_eachmice{1};
B=image_max_eachmice{1};
A(A<threshold)=0;
B(B<threshold)=0;
p_map_v=ds.image_diff(A,B,1,1);

A=image_max_base_eachmice{2};
B=image_max_eachmice{2};
A(A<threshold)=0;
B(B<threshold)=0;
p_map_a=ds.image_diff(A,B,1,1);

% A=image_max_eachmice{1};
% B=image_max_eachmice{2};
% A(A<threshold)=0;
% B(B<threshold)=0;
% p_map_AV=ds.image_diff(A,B,0,1);
% 
% figure;
% imagesc(double(p_map_AV>0.95)-1*double(p_map_AV<0.05))
% % imagesc(p_map_AV<0.05)
% 
% ap.wf_draw('ccf', [0.5 0.5 0.5]);
% axis image off

empty_area=zeros(450,426);
empty_area(double(p_map_v>0.95)+double(p_map_a>0.95)==2)=3;
empty_area(double(p_map_v>0.95)-double(p_map_a>0.95)==1)=2;
empty_area(-double(p_map_v>0.95)+double(p_map_a>0.95)==1)=1;


% drawnow
figure('Position',[50 50 400 300])
imagesc(empty_area)
my_colormap1 = [1 1 1;    % 0 = white
    1 0.5 0.5;    % 1 = blue
    0.5 0.5 1;    % 2 = red
    0.3 0.3 0.3];
colormap(my_colormap1);
ap.wf_draw('ccf', [0.5 0.5 0.5]);
axis image off

for curr_roi=[1 3]
B = bwboundaries(roi1(curr_roi).data.mask);             % 获取边界点
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2); % 画红色边线
end
 text(0.2, 0.8,roi1(1).name , ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',10, ...
                'FontWeight','normal');
 text(0.25, 0.9,roi1(3).name , ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',10, ...
                'FontWeight','normal');

clim([0 3])
ap.wf_draw('ccf', [0.5 0.5 0.5]);
cb=colorbar('eastoutside','Ticks', [0.375 1.125 1.875 2.625], ...
    'TickLabels', {'none','auditory','visual','both'});  % 可自定义标签
pos = cb.Position;      % [x, y, width, height]
pos(1) = 0.85;           % 右移
pos(2) = 0.22;
pos(3) = 0.03;           % 右移
pos(4) = 0.4;           % 变窄
cb.Position = pos;       % 应用修改

drawnow;
ap.prettyfig

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 1e.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});

%% fig 1f  task kernels traces
main_preload_vars = who;
load_dataset='wf_task_kernels';
load(fullfile(Path,'data',load_dataset));

tem_image_3=cellfun(@(x) permute(nanmean(plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[7 8],:)),4),[1,2,3,5,4]),...
    wf_task_kernels_across_day,'UniformOutput',false);

temp_each_roi=cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),tem_image_3,'UniformOutput',false);

buf3_roi_mean=cellfun(@(x)   nanmean(x,3) ,temp_each_roi,'UniformOutput',false )
buf3_roi_error=cellfun(@(x)  std(x,0,3,"omitmissing")./sqrt(size(x,3)) ,temp_each_roi,'UniformOutput',false  );

scale_mpfc=cellfun(@(x,y) [min(x([1:2],:),[],'all')-max(y([1:2],:),[],'all') max(x([1:2],:),[],'all')+max(y([1:2],:),[],'all');...
    min(x([1:2],:),[],'all')-max(y([1:2],:),[],'all') max(x([1:2],:),[],'all')+max(y([1:2],:),[],'all')],...
    buf3_roi_mean,buf3_roi_error,'UniformOutput',false);

scale_mpfc=[scale_mpfc(1);scale_mpfc(1)];
scale_apfc=cellfun(@(x,y)[min(x([3:4],:),[],'all')-max(y([3:4],:),[],'all') max(x([3:4],:),[],'all')+max(y([3:4],:),[],'all');...
    min(x([3:4],:),[],'all')-max(y([3:4],:),[],'all') max(x([3:4],:),[],'all')+max(y([3:4],:),[],'all')],...
    buf3_roi_mean,buf3_roi_error,'UniformOutput',false);
scale_sensory{1}= [min(buf3_roi_mean{1}(11,:),[],'all')-max(buf3_roi_error{1}(11,:),[],'all')...
    max(buf3_roi_mean{1}(11,:),[],'all')+max(buf3_roi_error{1}(11,:),[],'all')];
scale_sensory{2}= [min(buf3_roi_mean{2}(9,:),[],'all')-max(buf3_roi_error{2}(9,:),[],'all') ...
    max(buf3_roi_mean{2}(9,:),[],'all')+max(buf3_roi_error{2}(9,:),[],'all')];
scale_all=cellfun(@(x,y,z)  [x ;y ;z],  scale_sensory', scale_mpfc ,scale_apfc,'UniformOutput',false   );


[~,firing_begin_time]=cellfun(@(x) max(diff (x,1,2),[],2)   ,temp_each_roi ,'UniformOutput',false);
firng_begin_mean=cellfun(@(x) arrayfun(@(id)  nanmean(t_kernels(squeeze(x(id,:,:)))) ,1:length(roi1),'UniformOutput',true),...
    firing_begin_time,'UniformOutput',false );
firng_begin_error=cellfun(@(x) arrayfun(@(id)  std(t_kernels(squeeze(x(id,:,:))))/sqrt(size(x,3)) ,1:length(roi1),'UniformOutput',true),...
    firing_begin_time,'UniformOutput',false );



% drawnow

face_colors={[0 0 1],[1 0 0]};
select_area={[11  1 2 3 4 ];[9  1 2 3 4]};
figure('Position',[50 50 400 300])
mainfig=tiledlayout( 1,2, ...
    'TileSpacing', 'tight', 'Padding', 'none');
plot_fig=tiledlayout(mainfig,length(select_area{1}), 1, ...
    'TileSpacing', 'none', 'Padding', 'none');
plot_fig.Layout.Tile = 1;  % 明确放在主 layout 的第 1 个 tile
% plot_fig=cell(2,1);
for curr_group=1:2
    % plot_fig{curr_group}=tiledlayout(mainfig,length(select_area{1}), 1, ...
    %     'TileSpacing', 'none', 'Padding', 'none');
    % plot_fig{curr_group}.Layout.Tile = curr_group;  % 明确放在主 layout 的第 1 个 tile

    for curr_area =1:length(select_area{curr_group})
        a4=nexttile(plot_fig,curr_area)
        switch curr_group
            case 1
                yyaxis left
            case 2
                yyaxis right
        end
        hold on
        ap.errorfill(t_kernels,buf3_roi_mean{curr_group}(select_area{curr_group}(curr_area),:),...
            buf3_roi_error{curr_group}(select_area{curr_group}(curr_area),:),...
            face_colors{curr_group},0.1,1,2)

        plot(t_kernels,buf3_roi_mean{curr_group}(select_area{curr_group}(curr_area),:),...
            'Color',face_colors{curr_group},'LineWidth',2)
        xlim([-0.05 0.3])
        ylim([scale_all{curr_group}(curr_area,1) scale_all{curr_group}(curr_area,2)] )
        % ylim([-1.2 3.5]*1e-4)
        xline(0)
        axis off
        if curr_group==2
            if curr_area==1
                temp_name='sensory-L';
            else
                temp_name=roi1(select_area{curr_group}(curr_area)).name;
            end
            text(0.25, 0.9,temp_name , ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',10, ...
                'FontWeight','normal', 'Interpreter', 'none');

        end

    end
end
        aline=nexttile(plot_fig,length(select_area{curr_group}))
        line_loc=scale_all{curr_group}(curr_area,1)
        line([-0.05 0],[line_loc line_loc],'Color',[0 0 0],'LineStyle','-')
        
        line([-0.05 -0.05],[line_loc line_loc+1e-4],'Color',[0 0 0],'LineStyle','-')
            
        text(0.15, -0.1,'0.05s' ,  'Units','normalized', 'HorizontalAlignment','right', ...
                'VerticalAlignment','top', 'FontSize',8, 'FontWeight','normal');
     text(-0.15, 0.7,'10^{-4}\DeltaF/F_{0}' ,  'Units','normalized', 'HorizontalAlignment','right', ...
                'VerticalAlignment','top', 'FontSize',8, 'FontWeight','normal', 'Rotation',90);

nexttile(mainfig,2)
seq = [11  9  1 2 3 4];
hold on
colors = {[0 0 1], [1 0 0]};
group_defs = {
    [1 3 5 6], [1 2 4 5];  % group 1: use_seq, use_seq1
    [2 5 6],   [ 1 4 5];     % group 2: use_seq, use_seq1
    };
for curr_group = 1:2
    use_seq = group_defs{curr_group, 1};
    use_seq1 = group_defs{curr_group, 2};

    used_area = seq(use_seq);
    y_vals = use_seq1;
    x_vals = firng_begin_mean{curr_group}(used_area);
    x_errs = firng_begin_error{curr_group}(used_area);

    % 水平误差线
    line([x_vals - x_errs; x_vals + x_errs], [y_vals; y_vals], ...
        'Color', colors{curr_group}, 'LineWidth', 2)

    % 连线
    plot(x_vals, y_vals, '-o', 'Color', colors{curr_group}, 'LineWidth', 2,'MarkerSize',3,'MarkerFaceColor', colors{curr_group})

    % 散点
      % scatter(x_vals, y_vals, 20, colors{curr_group}, 'filled');


end

ylim([0.5 5.5])
yticks(1:6)
yticklabels({})

% yticklabels({'l-sensory', roi1(1).name, roi1(2).name, roi1(3).name, roi1(4).name})
xlabel('firing initial time (s)')
xlim([0.04 0.14])
xticks([0.04 0.09 0.14])
set(gca, 'YDir', 'reverse','FontSize',7.5,'Color','none')
box off
set(gcf,'Name', 'figure 1f')
set(gca,'YDir','reverse','Color','none')

temp_dis= cellfun(@(x,y)   permute(x,[1,3,2])   ,firing_begin_time,'UniformOutput',false  )
p_test=arrayfun(@(idx) ds.shuffle_test(temp_dis{1}(idx,:),temp_dis{2}(idx,:),0,1) ,3:4,'UniformOutput',true)
  ap.prettyfig;
%
exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 1f.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
clearvars('-except',main_preload_vars{:});

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
% clearvars('-except',main_preload_vars{:});

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

    ds.make_bar_plot(temp_plot_data([1 3 2 4]),colors1,'BarAlpha',1,'ShowDots',0,'CentralTendency','median','ShowErrorCaps',0)
        ds.make_bar_plot(temp_plot_data([1 3 2 4]),colors1,'BarAlpha',1,'ShowDots',0,'CentralTendency','median','ShowErrorCaps',0)

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


% fig = gcf;
% set(findall(fig, '-property', 'FontSize'), 'FontSize', get(gca,'FontSize') + 3);


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 3.eps'), ...
    'ContentType','vector');          % 导出为 EPS 矢量
% clearvars('-except',main_preload_vars{:});


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
% clearvars('-except',main_preload_vars{:});


%% Fig 5 passive kernels cross modality
main_preload_vars = who;
load(fullfile(Path,'data','wf_passive_kernels.mat'));
% kernels
tem_image=cellfun(@(x) cellfun(@(a) plab.wf.svd2px(U_master(:,:,1:size(a,1)),a(:,:,:,[ 10 11 16 17],:)),...
    x,'UniformOutput',false), wf_passive_kernels_across_day,'UniformOutput',false);
image_max_eachmice=cellfun(@(x,id) cellfun(@(a)    cat(3,  nanmean(max(a(:,:,kernels_period,id,[1 2],:),[],3),[5 6]),...
    nanmean(max(a(:,:,kernels_period,id,[3 4],:),[],3),[5 6])), x,'UniformOutput',false),tem_image,{3;2},'UniformOutput',false  );

tem_image_across_day=cellfun(@(x) cellfun(@(a) plab.wf.svd2px(U_master(:,:,1:size(a,1)),a(:,:,:,[4:11 18:23],:)),...
    x,'UniformOutput',false), wf_passive_kernels_across_day,'UniformOutput',false);
tem_roi_across_day= cellfun(@(x) cellfun(@(a) ds.make_each_roi(permute(max(a(:,:,kernels_period,:,:,:),[],3),[1,2,5,4,6,3]) ,size(a,5),roi1),...
    x,'UniformOutput',false),tem_image_across_day,'UniformOutput',false)
tem_roi_across_day_2=cellfun(@(x,id)   cellfun(@(a)    cat(3, a(:,:,id(1),:), nanmean(a(:,:,[id(2:3)],:),3) ),...
    x,'UniformOutput',false)  ,tem_roi_across_day,{[3,1,2];[2,1,3]},'UniformOutput',false)
tem_roi_across_day_mean=cellfun(@(x) cellfun(@(a) nanmean(a,4) , ...
    x,'UniformOutput',false),tem_roi_across_day_2,'UniformOutput',false );
tem_roi_across_day_error=cellfun(@(x) cellfun(@(a) std(a,0,4,'omitmissing')./sqrt(size(a,4)),...
    x,'UniformOutput',false) ,tem_roi_across_day_2,'UniformOutput',false );

scale=0.0002;
Color={'G','P'};

figure('Position', [2250 150 420 200]);
tiledlayout(2,4,'TileSpacing','tight', 'Padding', 'tight')
a1=cell(2,4);
for curr_passive=1:2
    for curr_group=1:2
        for curr_mod=1:2
            a1{curr_passive,2*curr_group-2+curr_mod}=nexttile
            imagesc( image_max_eachmice{curr_passive}{curr_group}(:,:,curr_mod))
            axis image off;
            colormap(a1{curr_passive,2*curr_group-2+curr_mod},ap.colormap(['W' Color{curr_group}] ));
            ap.wf_draw('ccf', [0.5 0.5 0.5]);
            clim(scale .* [0, 1]);
        end
    end
end

for curr_group=1:2
    cb=colorbar(a1{2,2*curr_group},'southoutside');
    % cb.Units = 'normalized';
    temp_pos= cb.Position;
    cb.Position = [ temp_pos(1)+0.02,0.08,0.1, 0.03];
    % cb.Position
end



exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 5_1.eps'), ...
    'ContentType','vector');
colors={  [  84 130 53 ]./255 ,[0.5 0.5 0.5];...
    [112  48 160]./255,[0.5 0.5 0.5]   }
%
figure('Position', [2250 250 440 220]);
plot_fig=tiledlayout(2,4,'TileSpacing','tight', 'Padding', 'tight')
vals=cell(2,1);
for curr_passive=1:2
    for curr_group=1:2
        temp_roi=0
        for use_roi=[1 3]
            temp_roi=temp_roi+1
            ax{curr_passive,2*curr_group-2+temp_roi}= nexttile(plot_fig,curr_passive*4-4+ 2*curr_group-2+temp_roi);
            hold on
            curr_statis=permute(tem_roi_across_day_2{curr_passive}{curr_group}(use_roi,:,:,:),[3,4,2,1]);


            for curr_stim=1:2
                color = colors{curr_group, curr_stim};
                ap.errorfill(1:5,tem_roi_across_day_mean{curr_passive}{curr_group}(use_roi,4:8,curr_stim),...
                    tem_roi_across_day_error{curr_passive}{curr_group}(use_roi,4:8,curr_stim),color,0.1 );

                ap.errorfill(6:10,tem_roi_across_day_mean{curr_passive}{curr_group}(use_roi,9:13,curr_stim),...
                    tem_roi_across_day_error{curr_passive}{curr_group}(use_roi,9:13,curr_stim),color,0.1 );
            end
            
            xline(5.5,'LineStyle',':')
axx = gca;
axx.YAxis.Exponent = -4;
            ylim(scale*[-0.1 1])
            xlim([0.5 10.5])
            xticks([1 5 6 10])
            if  use_roi==1
            ylabel( '\Delta F/F_{0}')
            end
            if curr_passive==2
                xticklabels({'-5','-1','0','4'})
                xlabel('day from transfer')
            else
                xticklabels([])

            end
            set(gca,'Color','none')
            temp{1}= nanmean(curr_statis(:,:,4:6),3);
            temp{2}= nanmean(curr_statis(:,:,7:8),3);
            temp{3}= nanmean(curr_statis(:,:,9:10),3);
            temp{4}= nanmean(curr_statis(:,:,11:13),3);

            vals{curr_group}{curr_passive}{temp_roi}=cellfun(@(x) ds.shuffle_test(x(1,:),x(2,:),0),temp,'UniformOutput',true)
                        temp_vals=cellfun(@(x) ds.shuffle_test(x(1,:),x(2,:),0),temp,'UniformOutput',true)

            
            xStart = [1 4 6 8 ]; xEnd = [3 5 7 10 ];

            line_thres=[-0.00001 -0.00001];
            arrayfun(@(a,b) line([a b],line_thres,'Color','k') ,xStart(temp_vals<0.05), xEnd(temp_vals<0.05));
            arrayfun(@(i) text(mean([xStart(i) xEnd(i)]), 0, '*','Color','r', 'FontSize',12,'HorizontalAlignment','center'), find(temp_vals<0.05));

            if curr_passive==1 
h = findobj(gca, 'Type', 'Line');
% 如果你想根据“绘制顺序”选择第 2 和第 6 条：
h = flipud(h);  % 翻转顺序，让 h(1) 对应第一条画的线
legend_name={'task', 'non-task'}
legend(h([temp_roi*2]), legend_name{temp_roi},'Location','northoutside','box','off','Orientation', 'horizontal');
            end
        end
    end
end

for curr_passive=1:2
    for curr_roi=1:2
        if curr_roi==1
            use_roi=1;
        else
            use_roi=3;
        end
        pos_ax  = get(ax{1,curr_passive*2-2+curr_roi}, 'Position')  % [left bottom width height]
        % pos_parent = plot_fig.Position;   % normalized relative to figure
        % pos_ax_fig = [ pos_parent(1) + pos_ax(1)*pos_parent(3), ...
        %     pos_parent(2) + pos_ax(2)*pos_parent(4), ...
        %     pos_ax(3)*pos_parent(3), ...
        %     pos_ax(4)*pos_parent(4) ];

        % 计算 inset 的位置（嵌在当前 tile 的左上角）
        inset_width = 0.3* pos_ax(3);    % inset 占 tile 宽度的 30%
        inset_height = 0.3 * pos_ax(4);   % inset 占 tile 高度的 30%
        inset_left = pos_ax(1) - 0* pos_ax(3);  % tile 左侧偏右一点
        inset_bottom = pos_ax(2) + 0.75 * pos_ax(4); % tile 底部偏上
        insetAx = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
        imagesc(roi1(use_roi).data.mask )
        ap.wf_draw('ccf', [0.5 0.5 0.5]);
        axis image off
        ylim([0 200])
        xlim([20 220])
        clim( [ 0, 1]);
        colormap( insetAx,ap.colormap('WK'));
        uistack(insetAx, 'bottom');
    end
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 5_2.eps'), ...
    'ContentType','vector');
% clearvars('-except',main_preload_vars{:});


%% fig s1a
main_preload_vars = who;

load(fullfile(Path,'data','example_trace.mat'));

time_period=[  min(find(timelite.timestamps-(photodiode_on_times(20)+-0.2)>0)),...
    min(find(timelite.timestamps-(photodiode_off_times(20)+0.5)>0))];

reward_timeline =reward_thresh(time_period(1):time_period(2));  % 示例数据
% 找出所有为 1 的索引
idx_ones = find(reward_timeline == 1);
% 计算相邻 1 之间的间隔
gap_lengths = diff(idx_ones) - 1;
% 找出 gap < 20 的区间索引
valid = find(gap_lengths < 20 & gap_lengths > 0);
% 创建逻辑掩码
mask = false(size(reward_timeline));
% 用 arrayfun 给 mask 中对应位置赋值为 true
idx_ranges = arrayfun(@(i) idx_ones(i)+1:idx_ones(i+1)-1, valid, 'UniformOutput', false);
mask(cell2mat(idx_ranges')) = true;
% 应用掩码修改 vec
reward_timeline(mask) = 1;
%
line_width=1;
font_size=8;
figure('Position',[50 50 200 150]);
% t1 = tiledlayout(4, 1, 'TileSpacing', 'loose', 'Padding', 'loose');

hold on
plot( photodiode_trace(time_period(1):time_period(2))>3,'LineWidth',line_width,'Color','k')
plot( reward_timeline-1.1,'LineWidth',line_width,'Color','k')
% plot(wheel_move(time_period(1):time_period(2))-2.2,'LineWidth',line_width,'Color','k')
% ylim([-0.1 1.1])
hold on
wheel_vel=wheel_velocity(time_period(1):time_period(2));
wheel_vel_norm = (wheel_vel - min(wheel_vel)) / (max(wheel_vel) - min(wheel_vel));

plot(wheel_vel_norm-2.5,'LineWidth',line_width,'Color','k')

ylim([-3 1.5])
axis off
% plot(lick_thresh(time_period(1):time_period(2))-4.4,'LineWidth',line_width,'Color','k')

xline(find(photodiode_trace(time_period(1):time_period(2))>3,1,'first'),...
    'LineStyle','--','LineWidth',1, 'Color',[0.5 0.5 0.5])
xline(find(reward_timeline==1,1,'first'),'LineStyle','--','LineWidth',1, 'Color',[0.5 0.5 0.5])
xline(find(wheel_move(time_period(1):time_period(2))==1,1,'first'),...
    'LineStyle','--','LineWidth',1, 'Color',[0.5 0.5 0.5])

line([find(photodiode_trace(time_period(1):time_period(2))>3,1,'first'),...
    find(wheel_move(time_period(1):time_period(2))==1,1,'first')],...
    [1.45,1.45], 'Color',[0.2 0.8 0.2],'LineWidth',1.5,'LineStyle','-')

line([find(photodiode_trace(time_period(1):time_period(2))>3,1,'first'),...
    find(reward_timeline==1,1,'first')],[1.25,1.25],...
    'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-')


labels = {'stim', 'reward',  'wheel velocity'};
y_positions = [0.3, -0.7, -2.2];

cellfun(@(label, y) text(160, y, label, ...
    'FontSize', font_size, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle'), ...
    labels, num2cell(y_positions));
xlim([-190 length(reward_timeline)])

text(250, 1.8, 'reaction time', ...
    'FontSize', font_size, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','Color',[0.2 0.8 0.2]), ...
    text(1000, 1.6, 'reward time', ...
    'FontSize', font_size, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','Color',[0.5 0.5 0.5]), ...
% 
exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s1a.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});
% exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s2e.eps'), ...
%     'ContentType','vector');

%% fig s1 b
main_preload_vars = who;

load(fullfile(Path,'data','example_behaviors.mat'));

animals={'DS010','DS015'};
titles={'mouse 1','mouse 2'}
figure('Position',[50 50 500 300])
mainlayout= tiledlayout(1,2)

for curr_animal =1:2
    switch curr_animal
        case 1
            colors_group=[0 0 1];

        case 2
            colors_group=[1 0 0];

    end
    animal=animals{curr_animal};
    raw_data_behavior=behavior. (animals{curr_animal});
    stage=1
    matches=unique(raw_data_behavior.workflow_name,'stable')
    learned_days=raw_data_behavior.rxn_l_mad_p(ismember(raw_data_behavior.workflow_name,{matches{stage}}),1)<0.01;


    sublayout= tiledlayout(mainlayout,2,1);
    sublayout.Layout.Tile=curr_animal



    sgtitle(sublayout,titles{curr_animal},'FontWeight','normal')
    for curr_state=1
        switch curr_state
            case 1
                temp_data=raw_data_behavior.stim2move_times(ismember(raw_data_behavior.workflow_name,{matches{stage}}),1);
                % ylabel_name='reaction time(s)';
                ylabel_name='RT (s)';

            case 2
                temp_data=raw_data_behavior.stim2lastmove_times(ismember(raw_data_behavior.workflow_name,{matches{stage}}),1);
                ylabel_name='reaction time(s)';

        end

        mergedVector = vertcat(temp_data{:});
        indexCells = cellfun(@(x, i) repmat(i, size(x)), temp_data, ...
            num2cell(1:numel(temp_data))', ...
            'UniformOutput', false);
        temp_idx=vertcat(indexCells{:});
        [~, firstIdx] = unique(temp_idx, 'stable');
        numCells = numel(temp_data);

        [unique_vals, ~, groupID] = unique(temp_idx, 'stable');
        mid_indices = splitapply(@(x) x(ceil(numel(x)/2)), (1:numel(temp_idx))', groupID);
        colors=zeros(numCells,3);
        colors(learned_days == 0, :) = repmat([0 0 0], sum(learned_days == 0), 1);
        colors(learned_days == 1, :) = repmat(colors_group, sum(learned_days == 1), 1);
        nexttile(sublayout)

        hold on
        for i = 1:numCells
            idx = (temp_idx == i);
            scatter(find(idx), mergedVector(idx), 10, colors(i,:), 'filled')
        end
        xline(firstIdx-0.5,':k')
        xlim([0 length(mergedVector)])
        ylim([0.05 20])
        ylabel(ylabel_name)

        xticks([])
        xlabel('days')
        if curr_state<4
            set(gca, 'YScale', 'log');
            yticks([1e-2 1e-1 1 10])
        else
            ylim([0 10])
            yticks([0 10])

        end
        drawnow
    end
    nexttile(sublayout)
    hold on
    yyaxis left
    set(gca, 'YColor', [0 0 0])
    ylabel('mad')
    temp_data=raw_data_behavior.stim2lastmove_mad(ismember(raw_data_behavior.workflow_name,{matches{stage}}),1);
    temp_data_null=raw_data_behavior.stim2lastmove_mad_null(ismember(raw_data_behavior.workflow_name,{matches{stage}}),1);
    plot(1:length(temp_data),temp_data,'LineStyle','-','Color',[0 0 0])
    plot(1:length(temp_data),temp_data_null,'LineStyle','--','Color',[0 0 0])
    set(gca, 'YScale', 'log');

    yyaxis right
    set(gca, 'YColor', colors_group)
    perform=(temp_data_null-temp_data)./(temp_data_null+temp_data)
    plot(1:length(perform),perform,'Color',colors_group)
    xlim([1 length(temp_data)])
    ylabel('perform')
    if sum(learned_days)>0
        xline(find(learned_days==1,1)-0.5,'LineStyle','--')
    end
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s1b.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});
%% fig s1c
main_preload_vars = who;
load(fullfile(Path,'data','behavior.mat'));

asso_day_mod1=cellfun(@(x)   sum(cellfun(@(a)  length(a)   ,x.reaction_time{:,2:3} ,'UniformOutput',true),2),behavior_aligned,'UniformOutput',false);
reaction_time_mod1=cellfun(@(x) cellfun(@(a)   a(end) , x.reaction_time{:,4},'UniformOutput',true) ,behavior_aligned,'UniformOutput',false);
perform_mod1=cellfun(@(x) cellfun(@(a)   a(end) , x.performance{:,4},'UniformOutput',true) ,behavior_aligned,'UniformOutput',false);


y_label={'first assocation day','RT (s)','perfromance'}
barColors = [[   187 205 174]./255;[ 198 172 217]./255]; % 浅蓝、浅红
scatterColors = [[84 130 53]./255; [112  48 160 ]./255]; % 深蓝、深红

barColors = [[0.8 0.8 1];[ 1 0.8 0.8]]; % 浅蓝、浅红
scatterColors = [[0.1 0.1 1]; [1 0.1 0.1]]; % 深蓝、深红
yscale={[1 9],[0 0.3 ],[0 1]};

figure('Position',[50 50 500 200]);

tiledlayout(1,3)
p_shuff=cell(3,1)
for curr_stage=1:3
    switch  curr_stage
        case 1
            temp_dat=asso_day_mod1;
        case 2
          temp_dat=  reaction_time_mod1;
        case 3
           temp_dat=  perform_mod1;
    end


% 计算均值和标准误差
means = cellfun(@mean, temp_dat);
stds = cellfun(@std, temp_dat);
nSamples = cellfun(@length, temp_dat);
sem = stds ./ sqrt(nSamples);  % 计算标准误 SEM
nexttile
% 创建柱状图，并确保 `bar` 只返回一个 `Bar` 对象数组
hold on;
% barHandle = bar(1:2, means, 0.5, 'FaceColor', 'flat','EdgeColor','none'); % 'FaceColor' 只能用于单个柱子时指定
plot(1:2,means,'k.','MarkerSize', 10)

% for i = 1:2
%     barHandle.CData(i,:) = barColors(i,:);
% 
% end
% x = barHandle.XData;

% 添加误差条
errorbar(1:2, means, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5); % 黑色误差条
jitterRange = 0.3;

for i = 1:2
    yvals = temp_dat{i};
    % 按 Y 值分组，避免相同值重叠
    [uniqueY, ~, idxGroup] = unique(yvals);
    xi = nan(size(yvals));
    for g = 1:numel(uniqueY)
        inds = find(idxGroup == g);
        nG = numel(inds);
        if nG > 1
            % 等间距分布 + 一点随机噪声
            baseJitter = linspace(-jitterRange/2, jitterRange/2, nG);
            noise = (rand(1, nG) - 0.5) * (jitterRange / nG);
            xi(inds) = i + baseJitter + noise;
        else
            xi(inds) = i + (rand - 0.5) * jitterRange;
        end
    end

    scatter(xi, yvals, 30, ...
        'MarkerFaceColor', scatterColors(i, :), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.7);
end

 ylim(yscale{curr_stage})
 xlim([0 3])
% 美化图像
xticks([1 2]);
xticklabels({'Visual', 'Auditory'});
xtickangle(45); % 将x轴标签旋转45度

ylabel(y_label{curr_stage});

yl=ylim;
yticks([yl(1) yl(2) ])
y_offset = (yl(2) - yl(1)) * 0.05;  % 横线高度偏移比例
% 横线和星号 y 位置
y_star = max([temp_dat{1}; temp_dat{2}]) + y_offset;
  p=  ranksum(temp_dat{1}, temp_dat{2});
 p_shuff{curr_stage} =ds.shuffle_test(temp_dat{1}, temp_dat{2},0,1)
% 判定星号数量
if p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'ns';  % 可选
     % stars = num2str(p);  % 可选
end
  % 添加横线和星号
plot([1 2], [y_star y_star], 'k-', 'LineWidth', 1.2);  % 横线
text(1.5, y_star + y_offset * 2, stars, ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'normal');
grid off;
hold off;
set(gca,'Color','none')
end


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s1c.eps'), ...
    'ContentType','vector');
% clearvars('-except',main_preload_vars{:});


%% fig s2a

surround_samplerate = 35;
surround_window = [-0.2,1];
t_task = surround_window(1):1/surround_samplerate:surround_window(2);
period_task=find(t_task>0&t_task<0.2);


main_preload_vars = who;
load(fullfile(Path,'data','wf_task_average_encoding.mat'));
 load(fullfile(Path,'data','wf_passive_average.mat'));


tem_average=cellfun(@(x)  cellfun(@(a)plab.wf.svd2px(U_master(:,:,1:size(a,1)),a),...
    x,'UniformOutput',false), wf_task_average_aligned,'UniformOutput',false);
tem_average=cellfun(@(x) cat(4,x{:}),tem_average,'UniformOutput',false);
image_average=cellfun(@(x) nanmean(max(x(:,:,period_task,:),[],3),4),tem_average, 'UniformOutput',false);

tem_encoding=cellfun(@(x)  cellfun(@(a)plab.wf.svd2px(U_master(:,:,1:size(a,1)),a),...
    x,'UniformOutput',false), wf_task_encoding_aligned,'UniformOutput',false);
tem_encoding=cellfun(@(x) cat(4,x{:}),tem_encoding,'UniformOutput',false);
image_encoding=cellfun(@(x) nanmean(max(x(:,:,kernels_period,:),[],3),4),tem_encoding, 'UniformOutput',false);

figure('Position',[50 50 400 200]);
tiledlayout(1,2,'TileSpacing','tight')
scale_image=0.02;
Color={'B','R'};

for curr_group=1:2
    ax=nexttile
    axs(curr_group) = ax;
    imagesc(image_average{curr_group})
    axis image off;
    clim(scale_image .* [0.25, 0.75]);
    colormap(ax, ap.colormap(['W' Color{curr_group}] ));
    clim(scale_image .* [0, 1]);
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end
drawnow
gap = 0.01;                          % 与子图的间距
cbh = 0.05;                          % colorbar 的高度（归一化坐标）
shrink = 1/3;                        % 宽度比例

for r = 1:2
    ax = axs(r);                     % 每个子图
    p = ax.Position;                 % [x y w h] 归一化
    cb = colorbar(ax,'southoutside');% 横向放在底下
    cb.Units = 'normalized';
    
    w = p(3)*shrink;                 % colorbar 宽度 = 子图宽度的 1/3
    x = p(1) + (p(3)-w)/2;           % 居中放置
    y = p(2) - gap - cbh;            % 放在子图底下，留出间距
    
    cb.Position = [x, y, w, cbh];    % [x y w h]
end


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s2a.eps'), ...
    'ContentType','vector');

figure('Position',[50 50 400 200]);
tiledlayout(1,2,'TileSpacing','tight')
scale_image=0.01;
Color={'B','R'};

for curr_group=1:2
    ax=nexttile
    axs(curr_group) = ax;
    imagesc(image_encoding{curr_group})
    axis image off;
    clim(scale_image .* [0.25, 0.75]);
    colormap(ax, ap.colormap(['W' Color{curr_group}] ));
    clim(scale_image .* [0, 1]);
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end
drawnow
gap = 0.01;                          % 与子图的间距
cbh = 0.05;                          % colorbar 的高度（归一化坐标）
shrink = 1/3;                        % 宽度比例

for r = 1:2
    ax = axs(r);                     % 每个子图
    p = ax.Position;                 % [x y w h] 归一化
    cb = colorbar(ax,'southoutside');% 横向放在底下
    cb.Units = 'normalized';
    
    w = p(3)*shrink;                 % colorbar 宽度 = 子图宽度的 1/3
    x = p(1) + (p(3)-w)/2;           % 居中放置
    y = p(2) - gap - cbh;            % 放在子图底下，留出间距
    
    cb.Position = [x, y, w, cbh];    % [x y w h]
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s2b.eps'), ...
    'ContentType','vector');


surround_samplerate = 35;
surround_window_passive = [-0.5,1];
t_passive = surround_window_passive   (1):1/surround_samplerate:surround_window_passive(2);
passive_boundary=0.2;
period_passive=find(t_passive>0&t_passive<passive_boundary); 

tem_passive_average=structfun(@(x)  cellfun(@(a)plab.wf.svd2px(U_master(:,:,1:size(a,1)),a),...
    x,'UniformOutput',false), wf_passive_average,'UniformOutput',false);
tem_passive_average=structfun(@(x) cat(5,x{:}),tem_passive_average,'UniformOutput',false);
image_passive_average=structfun(@(x) permute( nanmean(max(x(:,:,period_passive,:,:),[],3),5),[1,2,4,3,5]),tem_passive_average, 'UniformOutput',false);

scale_image=0.004;
Color={'B','R'};
figure('Position',[50 50 1200 200])
tiledlayout(1,6)
for curr_group=1:3
    ax=nexttile
    axs(curr_group) = ax;
    imagesc(image_passive_average.lcr_passive(:,:,curr_group))
    axis image off;
    colormap(ax, ap.colormap(['WB' ] ));
    clim(scale_image .* [0, 1]);
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end
for curr_group=1:3
    ax=nexttile
    axs(curr_group+3) = ax;
    imagesc(image_passive_average.hml_passive_audio(:,:,curr_group))
    axis image off;
    colormap(ax, ap.colormap(['WR' ] ));
    clim(scale_image .* [0, 1]);
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end

drawnow
gap = 0.01;                          % 与子图的间距
cbh = 0.01;                          % colorbar 的高度（归一化坐标）
shrink = 1/3;                        % 宽度比例
for r = [3 6]
    ax = axs(r);                     % 每个子图
    p = ax.Position;                 % [x y w h] 归一化
    cb = colorbar(ax,'eastoutside');% 横向放在底下
    cb.Units = 'normalized';
    
    w = p(4)*shrink;                 % colorbar 宽度 = 子图宽度的 1/3
    x = p(1)+0.12 ;           % 居中放置
    y = p(2) ;            % 放在子图底下，留出间距
    
    cb.Position = [x, y, cbh, w];    % [x y w h]
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s2c.eps'), ...
    'ContentType','vector');

clearvars('-except',main_preload_vars{:});


%% fig s3  kernels across time
main_preload_vars = who;

load(fullfile(Path,'data','wf_task_kernels'));
load(fullfile(Path,'data','wf_passive_kernels'));



tem_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,[7 8],:)),  wf_task_kernels_across_day,'UniformOutput',false);
image_acorss_time=cellfun(@(x)    nanmean(x,[4 5]),tem_image,'UniformOutput',false  );

tem_image_pass=cellfun(@(x,id) plab.wf.svd2px(U_master(:,:,1:size(x{id},1)),x{id}(:,:,:,[10 11],:)), ...
    wf_passive_kernels_across_day,{1;2},'UniformOutput',false);

image_pass_acorss_time=cellfun(@(x)    nanmean(x,[5 6]),tem_image_pass,'UniformOutput',false  );


scale_image=0.0003;
Color={'B','R'};

figure('Position', [50 50 900 700] )
mainfig=tiledlayout(8,1,'TileSpacing','none')
for curr_group=1:2
    subfig=tiledlayout(mainfig,1,sum(t_kernels>-0.05& t_kernels<0.25),'TileSpacing','none')
    subfig.Layout.Tile=curr_group;
    for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
        ax=nexttile(subfig)
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



% figure('Position', [50 50 900 200] )
% mainfig=tiledlayout(2,1,'TileSpacing','none')
for curr_group=1:2
    subfig=tiledlayout(mainfig,3,sum(t_kernels>-0.05& t_kernels<0.25),'TileSpacing','none')
    subfig.Layout.Tile=3*curr_group;
    subfig.Layout.TileSpan = [3 1];      % 跨 2 行 × 1 列

    for curr_stim=1:3
    for curr_frame=find(t_kernels>-0.05& t_kernels<0.25)
        ax=nexttile(subfig)
        imagesc(image_pass_acorss_time{curr_group}(:,:,curr_frame,curr_stim))
        axis image off;
        clim(scale_image .* [0, 1]);
        colormap(ax, ap.colormap(['W' Color{curr_group}] ));
        clim(scale_image .* [0, 1]);
        ap.wf_draw('ccf', [0.5 0.5 0.5]);
        % if curr_group==1&curr_stim==1
        % title( num2str(t_kernels(curr_frame),'%.2f'),'FontWeight','normal')
        % end
    end
    end
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s3.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});


%% Fig s4 a behaviors in cross modality task

% behavior
load(fullfile(Path,'data','behavior'));

reaction_time_1=cellfun(@(x) structfun(@(a) a([1:8 15:end],:),x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false)
reaction_time_mean=cellfun(@(x) structfun(@(a) nanmean(a([1:8 15:end],:),2),x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false)
reaction_time_error=cellfun(@(x) structfun(@(a) std(a([1:8 15:end],:),0,2,'omitmissing')./sqrt(size(a,2))  ,x,'UniformOutput',false) ,behavior_across_day,'UniformOutput',false)
perform_p= cellfun(@(x) structfun(@(a) ds.shuffle_test (nanmean(a(7:8,:),1),nanmean(a(9:10,:),1),1,2 )>0.95 ,x,'UniformOutput',false)  , reaction_time_1,'UniformOutput',false)
perform_pval= cellfun(@(x) structfun(@(a) 1-ds.shuffle_test (nanmean(a(7:8,:),1),nanmean(a(12:13,:),1),1,2 ) ,x,'UniformOutput',false)  , reaction_time_1,'UniformOutput',false)


colors = [ ...
    84 130 53  % #548235
    112  48 160  % #7030A0
    ] / 255;
figure('Position',[50 50 600 150])
tiledlayout(1,3,'TileSpacing','tight')
behav_para={'performance','itimove','velocity'}
behav_name={'Performance','Uncued/cued move','velocity'}
y_scale={[-0.1 1],[0 5],[0 3000]}
y_ticks={[0 1],[0 5],[0 3000]}
x_lim={[3 13],[3 13],[3 13]}
    style={'-','--'}

for curr_fig=1:3
nexttile
for curr_group=1:2
    hold on
    ap.errorfill(1:8, reaction_time_mean{curr_group}.(behav_para{curr_fig})(1:8),...
        reaction_time_error{curr_group}.(behav_para{curr_fig})(1:8),colors(curr_group,:),0.1,0)
    ap.errorfill(9:14, reaction_time_mean{curr_group}.(behav_para{curr_fig})(9:14),...
        reaction_time_error{curr_group}.(behav_para{curr_fig})(9:14),colors(curr_group,:),0.1,0)
    
   h(1)= plot(1:8, reaction_time_mean{curr_group}.(behav_para{curr_fig})(1:8)',...
        'Color',colors(curr_group,:),'LineStyle',style{curr_group},'LineWidth',2)
  h(2)= plot(9:14, reaction_time_mean{curr_group}.(behav_para{curr_fig})(9:14)',...
      'Color', colors(curr_group,:),'LineStyle',style{3-curr_group},'LineWidth',2)

    set(gca,'Color','none')
end

xline(8.5,'LineStyle','-','LineWidth',1,'Color',[0.5 0.5 0.5])
xlim(x_lim{curr_fig})
xticks([ 3 8 9 13 ])
xticklabels({'-6','-1','0','4'})
ylabel(behav_name{curr_fig})
ylim(y_scale{curr_fig})
yticks(y_ticks{curr_fig})
if curr_fig==3
yticklabels({'0','max'})
end
xlabel('day from transfer')
% test_p 是 [val1 val2]
offset = 0;  % 用于竖直堆叠
y_base = y_scale{curr_fig}(2);   % 星号基准高度，可以根据数据调节
test_p=cellfun(@(x)   x.(behav_para{curr_fig}),perform_p,'UniformOutput',true)
if test_p(1)
    text(9.5, y_base + offset, '*', 'Color',colors(1,:), 'FontSize',14, ...
        'HorizontalAlignment','center');
    offset = offset + y_scale{curr_fig}(2)*0.1;  % 往上移一格
end
if test_p(2)
    text(9.5, y_base + offset, '*', 'Color',colors(2,:), 'FontSize',14, ...
        'HorizontalAlignment','center');
end
set(gca,'Color','none')


end

% legend([h(1), h(2)], {'auditory', 'visual'} ,'Location','northwest','Box','off');


exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s4a.eps'), ...
    'ContentType','vector');          
% clearvars('-except',main_preload_vars{:});


%% Fig s4 b behavioral results in mixed tasks

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
% clearvars('-except',main_preload_vars{:});

%% fig s4c

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


%% fig s4d

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
% clearvars('-except',main_preload_vars{:});
