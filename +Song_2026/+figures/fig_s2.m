%% fig s2a
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



clearvars('-except',main_preload_vars{:});
