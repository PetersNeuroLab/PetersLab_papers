%% Fig s4 a behaviors in cross modality task
main_preload_vars = who;

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



exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s4a.eps'), ...
    'ContentType','vector');
 clearvars('-except',main_preload_vars{:});