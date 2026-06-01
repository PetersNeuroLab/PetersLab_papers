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
