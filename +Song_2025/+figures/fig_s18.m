%% EDF 18 probe position
main_preload_vars=who;


animals= { 'DS007','DS010','AP021','DS011','AP022',...
    'DS000','DS014','DS015','DS016'};

anterior_learned_idx_VA={[2 4],  [2 4],  2,   [2 4],  [2 4], [    ],   [   ],  [   ], [   ]};
anterior_learned_idx_AV={[ ],   [   ],  [],  [ ],    [ ],   [2 4 ],   [2 4],  [2 4], [2 4]};


groups={'VA','AV'}
animals_number{1}={'mouse 1','mouse 2','mouse 3','mouse 4','mouse 5'};
animals_number{2}={'mouse 1','mouse 2','mouse 3','mouse 4'};

for curr_group=1:2
    switch curr_group
        case 1
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
        case 2
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
    end
    temp_probe_position=cell(length(used_animals),1);

    for curr_animal=1:length(used_animals)

        animal=used_animals{curr_animal};
        temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);
        temp_probe_position{curr_animal}=temp_file_name.probe_positions(used_animals_idx{curr_animal},1);

    end
    probe_position_all{curr_group}=temp_probe_position;
    animals_group{curr_group}=used_animals;

end

fig_name={'VA','AV'}
for curr_group=1:2

    obj=ap.ccf_draw
    obj.draw_name('caudoputamen')
    obj.ccf_fig.Position=[50 50 800 200]

    cmap = jet(length(probe_position_all{curr_group}));

    for curr_animal=1:length(probe_position_all{curr_group})

        for curr_probe=1:length(probe_position_all{curr_group}{curr_animal})
            preload_vars = who;

            if isempty(curr_probe)

                continue
            end

            probe_line=probe_position_all{curr_group}{curr_animal}{curr_probe}';

            % Draw probes on coronal + saggital
            line(obj.ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',cmap(curr_animal,:));
            line(obj.ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',cmap(curr_animal,:));
            line(obj.ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',cmap(curr_animal,:));
            line(obj.ccf_axes(4),probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
                'linewidth',2,'color',cmap(curr_animal,:))
        end
    end


    colormap(cmap);          % 设置 colormap
    % cb=colorbar('Ticks', linspace(0, 1, length(animals_group{curr_group})), ...
    %     'TickLabels', animals_number{curr_group});  % 可自定义标签
        cb=colorbar('Ticks', linspace(0, 1, length(animals_group{curr_group})), ...
        'TickLabels', {});  % 可自定义标签
    % 获取当前 colorbar 的位置
    pos = cb.Position;     % pos = [x y width height]
    % 缩小宽度为原来的 50%
    pos(3) = pos(3) * 0.5;
    pos(4) = pos(4) * 0.5;
    cb.Position = pos;
    % sgtitle(groups{curr_group})
    % exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s6' fig_name{curr_group} '.eps']), ...
    % 'ContentType','vector');
    saveas(gcf, fullfile(Path,['figures\eps\Fig s6' fig_name{curr_group} '.png']));   % 保存为 PNG 文件

end



 clearvars('-except',main_preload_vars{:});