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
