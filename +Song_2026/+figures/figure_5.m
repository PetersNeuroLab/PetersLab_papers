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

figure('Position', [50 50 420 200]);
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
figure('Position', [50 50 440 220]);
plot_fig=tiledlayout(2,4,'TileSpacing','tight', 'Padding', 'tight')
vals=cell(2,1);
ax=cell(2,4);
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
clearvars('-except',main_preload_vars{:});
