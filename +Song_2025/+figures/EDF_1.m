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
 clearvars('-except',main_preload_vars{:});
