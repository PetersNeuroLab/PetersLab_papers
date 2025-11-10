clear all

Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';

animals= { 'DS007','DS010','AP021','DS011','AP022',...
    'DS000','DS014','DS015','DS016'};

anterior_learned_idx_VA={[2 4],  [2 4],  2,   [2 4],  [2 4], [    ],   [   ],  [   ], [   ]};
anterior_learned_idx_AV={[ ],   [   ],  [],  [ ],    [ ],   [2 4 ],   [2 4],  [2 4], [2 4]};

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

baseline_t_stim = [-0.2,0];
response_t_stim = [0,0.2];
psth_use_t_stim = t_bins >= response_t_stim(1) & t_bins <= response_t_stim(2);

baseline_t_move = [-0.3,-0.1];
response_t_move = [-0.1,0.1];
psth_use_t_move = t_bins >= response_t_stim(1) & t_bins <= response_t_stim(2);

titles={'L','M','R passive','4k','8k passive','12k','R task','8k task','R task move','8k task move','move'};



%% paradigm
main_preload_vars=who;

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
obj.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
obj.st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

figure('Position',[50 50 200 200]);
hold on
for curr_view = 2
    curr_outline = bwboundaries(squeeze((max(obj.av,[],curr_view)) > 1));
    % (only plot largest outline)
    [~,curr_outline_idx] = max(cellfun(@length,curr_outline));
    curr_outline_reduced = reducepoly(curr_outline{curr_outline_idx});
    plot(  curr_outline_reduced(:,2), ...
        curr_outline_reduced(:,1),'k','linewidth',2);
    set(gca,'YDir','reverse')
    axis(gca,'equal','off')

end

structure_name={'primary visual area','primary auditory area','caudoputamen'};
plot_structure_color = {[0.5 0.5 1],[1 0.5 0.5],[0 0 0]};
for curr_area=1:3
    plot_structure = find(strcmpi(obj.st.safe_name,structure_name{curr_area}));
    plot_structure_id = obj.st.structure_id_path{plot_structure};
    plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
        obj.st.structure_id_path));

    % Get structure color and volume
    slice_spacing = 5;
    plot_ccf_volume = ismember(obj.av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);
    % Plot 2D structure
    for curr_view = 2
        curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
        cellfun(@(x) plot(x(:,2)*slice_spacing, ...
            x(:,1)*slice_spacing,'color',plot_structure_color{curr_area},'linewidth',2),curr_outline(1))

        cellfun(@(x) fill(x(:,2)*slice_spacing, ...
            x(:,1)*slice_spacing, ...
            plot_structure_color{curr_area}, ...  % 填充颜色
            'FaceAlpha',0.5, ...                 % 透明度 (0=全透明,1=不透明)
            'EdgeColor',plot_structure_color{curr_area}, ...
            'LineWidth',2), ...
            curr_outline(1));

    end
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6a_1.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});


%%  probes_corornal slices
main_preload_vars=who;
bregma=[520,44,570];
groups={'VA','AV'}
colors=[[84 130 53]./255;[112  48 160]./255];
probe_position_all=cell(2,1);

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
    probe_position_all{curr_group}=vertcat(temp_probe_position{:});

end

% porbes_range in AP
temp_pos=cell2mat(cellfun(@(x)  x(1,:), vertcat(probe_position_all{:}),'UniformOutput',false));
temp_mid=nanmean([min(temp_pos,[],'all') max(temp_pos,[],'all')]);
position_slices=[temp_mid-30 temp_mid+30];


allen_atlas_path = fileparts(which('template_volume_10um.npy'));
obj.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
obj.st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

structure_name='caudoputamen';
plot_structure = find(strcmpi(obj.st.safe_name,structure_name));
plot_structure_id = obj.st.structure_id_path{plot_structure};
plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
    obj.st.structure_id_path));

% Get structure color and volume
% structure_color = hex2dec(reshape(obj.st.color_hex_triplet{plot_structure},2,[])')./255;
structure_color=[0.5 0.5 0.5]
plot_ccf_volume = ismember(obj.av,plot_ccf_idx);

dist2bregma=(bregma(1)-nanmean(position_slices))./100 %% mm

figure('Position',[50 50 200 200]);
hold on
for curr_view = 1
    curr_outline_out = bwboundaries(squeeze((max(obj.av(position_slices(1):position_slices(2),:,1:end/2),[],curr_view)) > 1));
    % curr_outline_out = bwboundaries(squeeze((obj.av(nanmean(position_slices),:,1:end/2))> 1));

    % (only plot largest outline)
    [~,curr_outline_idx] = max(cellfun(@length,curr_outline_out));
    curr_outline_reduced = reducepoly(curr_outline_out{curr_outline_idx});
    plot( ...
        curr_outline_reduced(:,2), ...
        curr_outline_reduced(:,1),'k','linewidth',2);
    % (draw 1mm scalebar)
    % line([0,0],[0,100],'color','k','linewidth',2);
end
set(gca,'YDir','reverse')
axis(gca,'equal','off')
curr_outline_area = bwboundaries(squeeze(max(plot_ccf_volume(position_slices(1):position_slices(2),:,:),[],curr_view)));
plot( curr_outline_area{1}(:,2), curr_outline_area{1}(:,1), ...
    'Color', structure_color, 'LineWidth', 2);
title(['AP: ' num2str(dist2bregma, '%.2f') ' mm'],'FontWeight','normal')

for curr_group=1:2
    cellfun(@(x)  line(x(3,:),x(2,:),'linewidth',2,'color',colors(curr_group,:)) ,   probe_position_all{curr_group},'UniformOutput',false)
end

exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6a.eps'), ...
    'ContentType','vector');
clearvars('-except',main_preload_vars{:});

%% probe position
main_preload_vars=who;

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


%% psth of all cells in passive
main_preload_vars = who;

colors={[84 130 53]./255,[112  48 160]./255};
image_color={'G','P'};
p_val=0.95;
for state=1
    switch state
        case 1
            plot_stim=[3 5];
            max_num=750;
            yscale=[-0.1 1.2];
            bar_scale=[0 0.5];
            clim_value=[0,5];
        case 2
            plot_stim=[7  8];
            max_num=750;
            yscale=[-0.1 2.5];
            bar_scale=[0 1];
            clim_value=[0,5];
    end


    sorting_stim=[3 5];

    figure('Position',[50 50 300 600]);
    p_fraction=cell(2,1)
    for curr_stim=1:2

        curr_sorting=sorting_stim(curr_stim);
        curr_plot=plot_stim(curr_stim);

        proportion_response=cell(2,1);
        proportion_response_overlay=cell(2,1);

        for curr_group=1:2
            switch curr_group
                case 1
                    used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
                    used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
                case 2
                    used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
                    used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));

            end

            temp_single_idx=cell(length(used_animals),1);
            temp_single_plot=cell(length(used_animals),1);
            temp_probe_position=cell(length(used_animals),1);
            temp_response=cell(length(used_animals),1);
            temp_response_plot=cell(length(used_animals),1);
            temp_response_plot_2=cell(length(used_animals),1);

            %
            for curr_animal=1:length(used_animals)

                animal=used_animals{curr_animal};
                temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);

                temp_single_plot{curr_animal}=temp_file_name.plot_single(used_animals_idx{curr_animal},1);
                temp_single_idx{curr_animal}=temp_file_name.plot_idx(used_animals_idx{curr_animal},1);
                temp_probe_position{curr_animal}=temp_file_name.all_celltypes(used_animals_idx{curr_animal},1);
                temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
                temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);
                temp_response_plot_2{curr_animal}=temp_file_name.all_event_response_signle_neuron_h2(used_animals_idx{curr_animal},1);


            end


            used_cell_type=vertcat(temp_probe_position{:});
            response=vertcat(temp_response{:});


            used_response_plot=vertcat(temp_response_plot{:});
            used_filter_plot_1=cellfun(@(x,y,z)  x(z(:,curr_sorting)>p_val ,:,curr_plot)  ,...
                used_response_plot,used_cell_type,response,'UniformOutput',false);
            used_filter_response_1=cellfun(@(x,y,z)  z(z(:,curr_sorting)>p_val ,:)  ,...
                used_response_plot,used_cell_type,response,'UniformOutput',false);
            used_plot_all_selected_1=vertcat(used_filter_plot_1{:});


            used_response_plot_2=vertcat(temp_response_plot_2{:});
            used_filter_plot_2=cellfun(@(x,y,z)  x(z(:,curr_sorting)>p_val ,:,curr_plot)  ,...
                used_response_plot_2,used_cell_type,response,'UniformOutput',false);
            used_filter_response_2=cellfun(@(x,y,z)  z(z(:,curr_sorting)>p_val ,:)  ,...
                used_response_plot_2,used_cell_type,response,'UniformOutput',false);
            used_plot_all_selected_2=vertcat(used_filter_plot_2{:});



            used_filter_plot_all=cellfun(@(x,y,z)  x(: ,:,curr_plot)  ,...
                used_response_plot,used_cell_type,response,'UniformOutput',false);
            used_plot_all=vertcat(used_filter_plot_all{:});



            proportion_response{curr_group}= cellfun(@(x) size(x,1) , used_filter_plot_1,'UniformOutput',true)./...
                cellfun(@(x) size(x,1) , used_filter_plot_all,'UniformOutput',true);

            [~,max_idx]=max(used_plot_all_selected_1(: ,psth_use_t_stim),[],2);
            [~,sort_idx] = sortrows( max_idx,"ascend");


            ax=subplot(5,2,[4*curr_group+curr_stim-4 ,4*curr_group+curr_stim-2])
            % ax=subplot(4,3,[curr_group+6*curr_stim-6 ,curr_group+6*curr_stim-3])
            % imagesc(t_bins,[],smoothdata(used_plot_all_selected(sort_idx,:),1,'gaussian',20))
            imagesc(t_bins,[],used_plot_all_selected_1(sort_idx,:))

            % colorbar('southoutside')
            colormap(ax,ap.colormap(['W' image_color{curr_group}]));
            clim(clim_value);
            xlim([-0.1 0.5])
            xline(0,'LineStyle',':')
            currentAx = gca; % 获取当前轴
            subplotPosition = currentAx.Position; % 获取位置和大小
            maxh=subplotPosition(4);
            maxb=subplotPosition(2);

            subplotPosition(4)=maxh/max_num*size(used_plot_all_selected_1,1);
            subplotPosition(2)=maxb+maxh-maxh/max_num*size(used_plot_all_selected_1,1);
            ax.Position=subplotPosition;
            axis off
            hold on

            if curr_stim==2

                cb = colorbar('south');  % 横向放在下方
                cb.Position = [subplotPosition(1) subplotPosition(2)-0.05 0.2 0.01];
                cb.Label.String = '\DeltaFR/FR_{0}';   % 给 colorbar 加标签

            end
            if curr_stim==1 &curr_group==2

                pos = ax.Position;   % [x0 y0 width height]

                % 横线 (annotation 是 figure 坐标)
                annotation('line', [pos(1)-0.05 pos(1)], [pos(2)-0.02 pos(2)-0.02], ...
                    'Color','k','LineWidth',2);
                % 竖线
                annotation('line', [pos(1)-0.05 pos(1)-0.05], [pos(2)-0.02 pos(2)-0.02+100*maxh/max_num], ...
                    'Color','k','LineWidth',2);

                % 文字
                annotation('textbox', [pos(1)-0.1 pos(2)-0.07 0.2 0.05], ...
                    'String','0.5 s', 'EdgeColor','none', ...
                    'HorizontalAlignment','center');
                annotation('textbox', [pos(1)-0.05 pos(2) 0.3 0.05], ...
                    'String','100 neurons', 'EdgeColor','none', ...
                    'HorizontalAlignment','center','Rotation',90);

            end



        end
        ax=subplot(5,2,8+curr_stim)
        ds.make_bar_plot(proportion_response,colors,0.2,10);
        hold on

        p_fraction{curr_stim} =  ranksum(proportion_response{1}, proportion_response{2});

        if p_fraction{curr_stim} < 0.05
            stars = repmat('*',1,sum(p_fraction{curr_stim}<[0.05 0.01 0.001]));
            y_sig = max(vertcat(proportion_response{:})) + 0.05;
            plot(1:2, [1 1]*y_sig, 'k-');
            text(1.5, y_sig+0.02, stars, 'HorizontalAlignment','center');
        end
        xticklabels({})
        ylabel('fraction','FontWeight','normal')
        ylim(bar_scale)
        yticks(bar_scale)
        box off
        set(gca,'color','none')

        drawnow
    end


    switch state
        case 1
            exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6c_2.eps'), ...
                'ContentType','vector');
        case 2
            exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s7a.eps'), ...
                'ContentType','vector');
    end

end
 clearvars('-except',main_preload_vars{:});

%% group average in passive by depth
main_preload_vars = who;

groups={'VA','AV','VA_nA'}
% groups={'VA','VA_nA'}

titles={'L','M','V passive','4k','A passive','12k','V task','A task','iti move'};
p_val=0.95
for state=1
    switch state
        case 1
            all_stim=[3 5];
        case 2
            all_stim=[7 8];

    end

    % colors={[84 130 53]./255,[112  48 160]./255};
    colors={[0.3 0.3 1],[1 0.3 0.3]};

    colors1={[0.1706    0.1275    0.1165],[0.3294    0.5098    0.2078], [0.7451    0.8667    0.6706];...
        [0.2706    0.0353    0.4667],[0.4392    0.1882    0.6275], [0.8196    0.7216    0.9019]};
    colors_image={'G','P'}
    z_min = 0;
    z_max = 250;
    bin_size_z = 25; % 单位：μm，根据实际尺度调整

    z_edges = [z_min:bin_size_z:z_max,inf];

    fig2 = figure('Position',[50 50 400 400]);
    tl2 = tiledlayout(length(all_stim),3,'TileSpacing','loose');
    firing_rates_max_mice=cell(2,1);
    for curr_group=1:2

        switch curr_group
            case 1
                used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
                used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
            case 2
                used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
                used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
        end

        temp_single_idx=cell(length(used_animals),1);
        temp_single_plot=cell(length(used_animals),1);
        temp_probe_position=cell(length(used_animals),1);
        temp_response=cell(length(used_animals),1);
        temp_response_plot=cell(length(used_animals),1);
        temp_cell_position=cell(length(used_animals),1);
        for curr_animal=1:length(used_animals)
            animal=used_animals{curr_animal};
            temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);
            temp_single_plot{curr_animal}=temp_file_name.plot_single(used_animals_idx{curr_animal},1);
            temp_single_idx{curr_animal}=temp_file_name.plot_idx(used_animals_idx{curr_animal},1);
            temp_probe_position{curr_animal}=temp_file_name.all_celltypes(used_animals_idx{curr_animal},1);
            temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
            temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);
            temp_cell_position{curr_animal}= cellfun(@(x,y) x-y , temp_file_name.all_cell_ccf_position_sorted(used_animals_idx{curr_animal},1),...
                temp_file_name.striatal_surface_position(used_animals_idx{curr_animal},1),'UniformOutput',false);
        end

        single_neuron_each_rec_1=vertcat(temp_response_plot{:});
        single_neuron_all_plot=cat(1,single_neuron_each_rec_1{:});
        response_each_rec=vertcat(temp_response{:});
        response_all=cat(1,response_each_rec{:});
        single_neuron_each_position=vertcat(temp_cell_position{:});
        single_neuron_position_all=cat(1,single_neuron_each_position{:});

        neuron_count_map_all=cell(2,1);
        neuron_count_map=cell(2,1);
        neuron_count_map_overlay=cell(2,1);
        for curr_sorting=1:length(all_stim)



            used_stim=all_stim(curr_sorting);

            neuron_coords_all= cellfun(@(x)   x(:,2) ,single_neuron_each_position,'UniformOutput',false );
            neuron_coords_each= cellfun(@(x,y) x(y(:,used_stim)>p_val,2),single_neuron_each_position,response_each_rec,'UniformOutput',false)


            if used_stim==3|used_stim==5
                temp_single_idx=[3 5]
            else
                temp_single_idx=[7 8]
            end




            firing_rates=cellfun(@(x) x(:,:,used_stim) ,single_neuron_each_rec_1,'UniformOutput',false);

            % === Step 1: 投影到冠状面（y-z） ===
            projected_coords = neuron_coords_each; % 取 y 和 z
            projected_coords_all = neuron_coords_all; % 取 y 和 z


            % 分配每个神经元的 bin 索引
            [neuron_count_map_all{curr_sorting},~,binIdx_all] = cellfun(@(x) histcounts(x, z_edges),projected_coords_all,'UniformOutput',false);
            [neuron_count_map{curr_sorting},~,binIdx]= cellfun(@(x) histcounts(x, z_edges),projected_coords,'UniformOutput',false);



            firing_rates_bins =cellfun(@(x,idx) arrayfun(@(col) ...
                accumarray(idx, x(:,col), [length(z_edges)-1,1], @mean, NaN), ...
                1:size(x,2), 'UniformOutput', false),firing_rates,binIdx_all, 'UniformOutput', false);

            firing_rates_bins1=cellfun(@(x) cat(2,x{:}),firing_rates_bins,'UniformOutput',false);
            firing_rates_bins2=nanmean(cat(3,firing_rates_bins1{:}),3);


            figure(fig2)
            % ax2=nexttile(tl2,3*curr_sorting-3+curr_group)
            ax2=nexttile(tl2,curr_sorting-3+3*curr_group)

            h2=imagesc(t_bins,z_edges(1:end-1),firing_rates_bins2)
            xlim([-0.1 0.5]);
            hold on
            xline(0,'LineStyle',':')
            ylim([z_edges(1)-0.5*bin_size_z  z_edges(end-1)+0.5*bin_size_z])

            if used_stim<7
                clim([0 2])
            else
                clim([0 4])
            end
            colormap(ax2,ap.colormap(['W' colors_image{curr_group}]))
            if curr_group==2
                xlabel('time (s)');
                xticks([-0.1 0.5]);
            else
                xticks([]);

            end
            if curr_sorting==2
                cb = colorbar(ax2,'eastoutside');  % 横向放在下方
                pos = cb.Position;   % [left bottom width height]
                pos(4) = pos(4) - 0.2;   % 缩短高度
                pos(3) = pos(3) /2;   % 缩短高度

                cb.Position = pos;

                cb.Label.String = '\DeltaFR/FR_{0}';   % 给 colorbar 加标签

            end



            if curr_sorting==1
                % title(titles(used_stim),'FontWeight','normal')
                ylabel('depth (\mum)');
                yticks([z_edges(1) z_edges(end-1)]);
            else
                yticks([]);

            end


            % MUA max
            % ax1=nexttile(tl2,3*curr_sorting)
            ax1=nexttile(tl2,3*curr_group)

            firing_rates_max=cellfun(@(x) max(x(:,psth_use_t_stim),[],2),firing_rates_bins1,'UniformOutput',false );
            firing_rates_max_mice{curr_group}{curr_sorting}=cat(2,firing_rates_max{:});
            ap.errorfill(z_edges(1:end-1) , smoothdata(nanmean(cat(3,firing_rates_max{:}),3),'gaussian',4),...
                smoothdata( std(cat(3,firing_rates_max{:}),0,3,'omitmissing')./sqrt(size(cat(3,firing_rates_max{:}),3)),'gaussian',4),...
                colors{curr_sorting},0.1,0.5);

            if used_stim<7
                ylim(ax1,[0 3])
            else
                ylim(ax1,[0 8])
            end
            xlim(ax1,[z_edges(1) z_edges(end-1)])

            ylabel('\DeltaFR/FR_0')
          
            if curr_group==2
                line([75 125],[0.2 0.2],'color','r')
            end
            xticks([]);

            view(ax1,90, 90);
             set(gca, 'Color', 'none')  


        end


    end


    cellfun(@(x) arrayfun(@(id) ds.shuffle_test(x{1}(id,:),x{2}(id,:),1,1),...
        1:10,'UniformOutput',false),firing_rates_max_mice,'UniformOutput',false)




    switch state
        case 1
            exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6d_2.eps'), ...
                'ContentType','vector');
        case 2
            exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s7b.eps'), ...
                'ContentType','vector');
    end

end
 clearvars('-except',main_preload_vars{:});

%% scatter of visual vs auditory 
main_preload_vars = who;

colors={[84 130 53]./255,[112  48 160]./255};
image_color={'G','P'};
p_val=0.95;
peak=cell(2,1);
idx_temp=cell(2,1);
types_all=cell(2,1);

sorting_stim=[3 5];


for curr_group=1:2
    switch curr_group
        case 1
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
        case 2
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));

    end

    temp_single_idx=cell(length(used_animals),1);
    temp_single_plot=cell(length(used_animals),1);
    temp_probe_position=cell(length(used_animals),1);
    temp_response=cell(length(used_animals),1);
    temp_response_plot=cell(length(used_animals),1);
    temp_response_plot_2=cell(length(used_animals),1);
    temp_type=cell(length(used_animals),1);
    
    for curr_animal=1:length(used_animals)

        animal=used_animals{curr_animal};
        temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);

        temp_single_plot{curr_animal}=temp_file_name.plot_single(used_animals_idx{curr_animal},1);
        temp_single_idx{curr_animal}=temp_file_name.plot_idx(used_animals_idx{curr_animal},1);
        temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
        temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);
        temp_type{curr_animal}=temp_file_name.all_celltypes(used_animals_idx{curr_animal},1);

    end

    types=vertcat(temp_type{:});
    response=vertcat(temp_response{:});
    used_response_plot=vertcat(temp_response_plot{:});


    used_plot_all_selected_1=vertcat(used_response_plot{:});

    temp=reshape(zscore(reshape(used_plot_all_selected_1(: ,:,sorting_stim),size(used_plot_all_selected_1,1),[],1),0,1),...
        size(used_plot_all_selected_1,1),size(used_plot_all_selected_1,2),2);
    temp=used_plot_all_selected_1(: ,:,sorting_stim);
    [peak_v,max_idx_v]=max(temp(: ,psth_use_t_stim,1),[],2);
    [peak_a,max_idx_a]=max(temp(: ,psth_use_t_stim,2),[],2);
    peak{curr_group}={peak_v,peak_a};

    used_response_all=vertcat(response{:});
    temp_response= used_response_all(:,sorting_stim)>p_val;

    temp_data{1} = find(temp_response(:,1)==1 & temp_response(:,2)==0);
    temp_data{3} = find(temp_response(:,1)==1 & temp_response(:,2)==1);
    temp_data{2} = find(temp_response(:,1)==0 & temp_response(:,2)==1);
    idx_temp{curr_group}=temp_data;


    fields = {'msn','fsi','tan'};

  types_all{curr_group}=cellfun(@(name) feval(@(C) vertcat(C{:}),  arrayfun(@(id) types{id}.(name) ,1:length(types),'UniformOutput',false)),fields,'UniformOutput',false)

end



scales = {[0 20 80],[0 8 20],[0 4 10]};
colors = {[0 0 1],[1 0 0],[ 1.0, 0.647, 0.0]};
    title_name={'VA','AV'}

% for curr_type=1:3
% 
%     for curr_group=1:2
% 
% 
%         figure('Position',[50 50 200 200])
% 
%         ax1 = axes('Position',[0.2 0.2 0.5 0.5]); % 左边大一些
%         hold on
%         cellfun(@(id,color) scatter(peak{curr_group}{1}(intersect(find(types_all{curr_group}{curr_type}) , id)),...
%             peak{curr_group}{2}(intersect(find(types_all{curr_group}{curr_type}) , id)),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
%             idx_temp{curr_group},colors,'UniformOutput',false)
% 
% length(intersect(find(types_all{curr_group}{curr_type}) , idx_temp{curr_group}{3}))
% 
% 
%         xlim(scales{curr_type}(1:2))
%         ylim(scales{curr_type}(1:2))
%         xticks(scales{curr_type}(1:2))
%         yticks(scales{curr_type}(1:2))
%         set(gca,'Color','none')
%         xlabel('visual response (\DeltaFR/FR_{0})')
%         ylabel('auditory response (\DeltaFR/FR_{0})')
% 
%         ax2 = axes('Position',[0.75 0.2 0.2 0.5]);
%         hold on
% 
%         cellfun(@(id,color) scatter(peak{curr_group}{1}(intersect(find(types_all{curr_group}{curr_type}) , id)),...
%             peak{curr_group}{2}(intersect(find(types_all{curr_group}{curr_type}) , id)),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
%             idx_temp{curr_group},colors,'UniformOutput',false)
%         xlim(scales{curr_type}(2:3))
%         ylim(scales{curr_type}(1:2))
%         xticks(scales{curr_type}(3))
%         set(gca,'YTick',[])
%         set(gca,'YColor','none')
%         set(gca,'Color','none')
% 
%         ax3 = axes('Position',[0.2 0.75 0.5 0.2]);
%         hold on
% 
%         cellfun(@(id,color) scatter(peak{curr_group}{1}(intersect(find(types_all{curr_group}{curr_type}) , id)),...
%             peak{curr_group}{2}(intersect(find(types_all{curr_group}{curr_type}) , id)),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
%             idx_temp{curr_group},colors,'UniformOutput',false)
%         ylim(scales{curr_type}(2:3))
%         xlim(scales{curr_type}(1:2))
%         set(gca,'Color','none')
%         yticks(scales{curr_type}(3))
%         set(gca,'XTick',[])
%         set(gca,'XColor','none')
% 
%         ax4 = axes('Position',[0.75 0.75 0.2 0.2]);
%         hold on
% 
%         cellfun(@(id,color) scatter(peak{curr_group}{1}(intersect(find(types_all{curr_group}{curr_type}) , id)),...
%             peak{curr_group}{2}(intersect(find(types_all{curr_group}{curr_type}) , id)),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
%             idx_temp{curr_group},colors,'UniformOutput',false)
%         ylim(scales{curr_type}(2:3))
%         xlim(scales{curr_type}(2:3))
%         axis off
%         title(title_name{curr_group},'FontWeight','normal')
% 
%         annotation('line',[0.74 0.76],[0.19 0.21],'Color','k','LineWidth',1) % 左斜杠
%         annotation('line',[0.69 0.71],[0.19 0.21],'Color','k','LineWidth',1) % 右斜杠
%         annotation('line',[0.19 0.21],[0.74 0.76],'Color','k','LineWidth',1) % 左斜杠
%         annotation('line',[0.19 0.21],[0.69 0.71],'Color','k','LineWidth',1) % 右斜杠
% 
% 
% 
%     end
% 
% end
% 


scales_all=[0 20 80];

for curr_group=1:2
    figure('Position',[50 50 200 200])
    ax1 = axes('Position',[0.2 0.2 0.5 0.5]); % 左边大一些
    hold on
    cellfun(@(id,color) scatter(peak{curr_group}{1}(id),...
        peak{curr_group}{2}(id),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
        idx_temp{curr_group},colors,'UniformOutput',false)

    xlim(scales_all(1:2))
    ylim(scales_all(1:2))
    xticks(scales_all(1:2))
    yticks(scales_all(1:2))
    set(gca,'Color','none')
    xlabel({'visual response' ;'(\DeltaFR/FR_{0})'})
    ylabel({'auditory response';' (\DeltaFR/FR_{0})'})

    ax2 = axes('Position',[0.75 0.2 0.2 0.5]);
    hold on

    cellfun(@(id,color) scatter(peak{curr_group}{1}(id),...
        peak{curr_group}{2}(id),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
        idx_temp{curr_group},colors,'UniformOutput',false)
    xlim(scales_all(2:3))
    ylim(scales_all(1:2))
    xticks(scales_all(3))
    set(gca,'YTick',[])
    set(gca,'YColor','none')
    set(gca,'Color','none')

    ax3 = axes('Position',[0.2 0.75 0.5 0.2]);
    hold on
    cellfun(@(id,color) scatter(peak{curr_group}{1}(id),...
        peak{curr_group}{2}(id),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
        idx_temp{curr_group},colors,'UniformOutput',false)
    ylim(scales_all(2:3))
    xlim(scales_all(1:2))
    set(gca,'Color','none')
    yticks(scales_all(3))
    set(gca,'XTick',[])
    set(gca,'XColor','none')

    ax4 = axes('Position',[0.75 0.75 0.2 0.2]);
    hold on
    cellfun(@(id,color) scatter(peak{curr_group}{1}(id),...
        peak{curr_group}{2}(id),20,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5),...
        idx_temp{curr_group},colors,'UniformOutput',false)
    ylim(scales_all(2:3))
    xlim(scales_all(2:3))
    axis off
    title(title_name{curr_group},'FontWeight','normal')

    annotation('line',[0.74 0.76],[0.19 0.21],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.69 0.71],[0.19 0.21],'Color','k','LineWidth',1) % 右斜杠
    annotation('line',[0.19 0.21],[0.74 0.76],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.19 0.21],[0.69 0.71],'Color','k','LineWidth',1) % 右斜杠


exportgraphics(gcf, fullfile(Path,['figures\eps\Fig 6_scatter' num2str(curr_group)  '.eps']), ...
                'ContentType','vector');
end

 clearvars('-except',main_preload_vars{:});


%%  overlay
main_preload_vars = who;
p_val=0.95;
colors={[84 130 53]./255,[112  48 160]./255};

line_error=cell(2,1);
line_mean=cell(2,1);
response_both_fration=cell(2,1);


line_error_both=cell(2,1);
line_mean_both=cell(2,1);
response_both_fration_both=cell(2,1);

for curr_group=1:2
    switch curr_group
        case 1
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
        case 2
            used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
            used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
    end

    temp_response=cell(length(used_animals),1);
    temp_response_plot=cell(length(used_animals),1);

    for curr_animal=1:length(used_animals)
        animal=used_animals{curr_animal};
        temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);
        temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
        temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);

    end
    response= cellfun(@(x)  x(:,[3 5])>p_val,vertcat(temp_response{:}),'UniformOutput',false );
    ptsh= cellfun(@(x)  x(:,:,[3 5]),vertcat(temp_response_plot{:}),'UniformOutput',false );

    for curr_stage=1:2


        response_sort= cellfun(@(x)    sortrows( x( sum(x,2)>0&x(:,curr_stage)==1,:) ,[-1,-2] ) ,response,'UniformOutput',false );
        response_sort_all=cat(1,response_sort{:});
        response_both_fration{curr_group}{curr_stage}  = cellfun(@(x) sum(sum(x,2)==2)./size(x,1) ,...
            response_sort,'UniformOutput',true);
        response_overlap=sum(cellfun(@(x) sum(sum(x,2)==2) ,response_sort,'UniformOutput',true));
        n_shuff=1000;
        response_both_shuff=   arrayfun(@(shuff) ap.shake(response_sort_all,1),1:n_shuff,...
            'UniformOutput',false);
        response_both_fration_shuff=cellfun(@(x)   sum(x(1:response_overlap,3-curr_stage))/ length(x),...
            response_both_shuff,'UniformOutput',true);
        line_mean{curr_group}(curr_stage)= (prctile(response_both_fration_shuff,95)+ prctile(response_both_fration_shuff,5))/2;
        line_error{curr_group}(curr_stage)=(prctile(response_both_fration_shuff,95)- prctile(response_both_fration_shuff,5))/2;
    end

    response_sort_both= cellfun(@(x)    sortrows( x( sum(x,2)>0,:) ,[-1,-2] ) ,response,'UniformOutput',false );
    response_sort_all_both=cat(1,response_sort_both{:});
    % response_both_fration_v{curr_group}  = cellfun(@(x) sum(sum(x,2)==2)./size(x,1) ,...
    %     response_sort_both,'UniformOutput',true);
    response_both_fration_both{curr_group}  = cellfun(@(x) sum(sum(x,2)==2)./size(x,1) ,...
        response_sort_both,'UniformOutput',true);
    n_shuff=1000;
    response_both_shuff_both=   arrayfun(@(shuff) ap.shake(response_sort_all_both,1),1:n_shuff,...
        'UniformOutput',false);
    response_both_fration_shuff_both=cellfun(@(x)    sum( sum(x,2)==2)/ length(x),...
        response_both_shuff_both,'UniformOutput',true);
    line_mean_both{curr_group}= (prctile(response_both_fration_shuff_both,95)+ prctile(response_both_fration_shuff_both,5))/2;
    line_error_both{curr_group}=(prctile(response_both_fration_shuff_both,95)- prctile(response_both_fration_shuff_both,5))/2;



end


figure('Position',[50 50 400 120]);
tiledlayout(1,3)
maker_size=20;
temp_name={'V','A'}
for curr_stage=1:2
    nexttile
    bar_temp={response_both_fration{1}{curr_stage},response_both_fration{2}{curr_stage}};
    ds.make_bar_plot(bar_temp,colors',0.2,maker_size);
    ylim([0 0.8])
    yticks([0 0.8])
    set(gca,'color','none')
    xticks([])
    xlabel(['overlay in ' temp_name{curr_stage}])
end

nexttile;
hold on

ds.make_bar_plot(response_both_fration_both,colors',0.2,maker_size);
arrayfun(@(curr_group)  ap.errorfill(curr_group-0.4:0.8:curr_group+0.4,[line_mean_both{curr_group} line_mean_both{curr_group}],...
    [line_error_both{curr_group} line_error_both{curr_group}],[0.5 0.5 0.5],1),1:2);
ylim([0 0.5])
yticks([0 0.5])
xticks([])
xlabel(['overlay in V∪A' ])
set(gca,'color','none')
exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6e_2.eps'), ...
    'ContentType','vector');

 clearvars('-except',main_preload_vars{:});

%%

% %%  psth of all cells in task
% main_preload_vars = who;
% 
% colors={[84 130 53]./255,[112  48 160]./255};
% image_color={'G','P'};
% p_val=0.95;
% plot_stim=[7  8];
% max_num=750;
% yscale=[-0.1 2.5];
% bar_scale=[0 1];
% clim_value=[0,5];
% 
% sorting_stim=[3 5];
% figure('Position',[50 50 400 300]);
% for curr_stim=1:2
% 
%     curr_sorting=sorting_stim(curr_stim);
%     curr_plot=plot_stim(curr_stim);
% 
%     proportion_response=cell(2,1);
%     proportion_response_overlay=cell(2,1);
% 
%     for curr_group=1:2
%         switch curr_group
%             case 1
%                 used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
%                 used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
%             case 2
%                 used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
%                 used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
% 
%         end
% 
%         temp_single_idx=cell(length(used_animals),1);
%         temp_single_plot=cell(length(used_animals),1);
%         temp_probe_position=cell(length(used_animals),1);
%         temp_response=cell(length(used_animals),1);
%         temp_response_plot=cell(length(used_animals),1);
%         temp_response_plot_2=cell(length(used_animals),1);
% 
%         %
%         for curr_animal=1:length(used_animals)
% 
%             animal=used_animals{curr_animal};
%             temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);
% 
%             temp_single_plot{curr_animal}=temp_file_name.plot_single(used_animals_idx{curr_animal},1);
%             temp_single_idx{curr_animal}=temp_file_name.plot_idx(used_animals_idx{curr_animal},1);
% 
%             temp_probe_position{curr_animal}=temp_file_name.all_celltypes(used_animals_idx{curr_animal},1);
%             temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
%             temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);
%             temp_response_plot_2{curr_animal}=temp_file_name.all_event_response_signle_neuron_h2(used_animals_idx{curr_animal},1);
% 
% 
%         end
% 
%         % used_single_plot=vertcat(temp_single_plot{:});
%         % used_single_idx=vertcat(temp_idx{:});
% 
%         used_cell_type=vertcat(temp_probe_position{:});
%         response=vertcat(temp_response{:});
% 
% 
%         used_response_plot=vertcat(temp_response_plot{:});
%         used_filter_plot_1=cellfun(@(x,y,z)  x(z(:,curr_sorting)>p_val ,:,curr_plot)  ,...
%             used_response_plot,used_cell_type,response,'UniformOutput',false);
%         used_filter_response_1=cellfun(@(x,y,z)  z(z(:,curr_sorting)>p_val ,:)  ,...
%             used_response_plot,used_cell_type,response,'UniformOutput',false);
%         used_plot_all_selected_1=vertcat(used_filter_plot_1{:});
% 
% 
%         used_response_plot_2=vertcat(temp_response_plot_2{:});
%         used_filter_plot_2=cellfun(@(x,y,z)  x(z(:,curr_sorting)>p_val ,:,curr_plot)  ,...
%             used_response_plot_2,used_cell_type,response,'UniformOutput',false);
%         used_filter_response_2=cellfun(@(x,y,z)  z(z(:,curr_sorting)>p_val ,:)  ,...
%             used_response_plot_2,used_cell_type,response,'UniformOutput',false);
%         used_plot_all_selected_2=vertcat(used_filter_plot_2{:});
% 
% 
% 
%         used_filter_plot_all=cellfun(@(x,y,z)  x(: ,:,curr_plot)  ,...
%             used_response_plot,used_cell_type,response,'UniformOutput',false);
%         used_plot_all=vertcat(used_filter_plot_all{:});
% 
% 
% 
%         proportion_response{curr_group}= cellfun(@(x) size(x,1) , used_filter_plot_1,'UniformOutput',true)./...
%             cellfun(@(x) size(x,1) , used_filter_plot_all,'UniformOutput',true);
% 
%         [~,max_idx]=max(used_plot_all_selected_1(: ,psth_use_t_stim),[],2);
%         [~,sort_idx] = sortrows( max_idx,"ascend");
% 
% 
%         % ax=subplot(6,2,[4*curr_group+curr_stim-4 ,4*curr_group+curr_stim-2])
%         ax=subplot(4,3,[curr_group+6*curr_stim-6 ,curr_group+6*curr_stim-3])
% 
%         % imagesc(t_bins,[],smoothdata(used_plot_all_selected(sort_idx,:),1,'gaussian',20))
%         imagesc(t_bins,[],used_plot_all_selected_1(sort_idx,:))
% 
%         % colorbar('southoutside')
%         colormap(ax,ap.colormap(['W' image_color{curr_group}]));
%         clim(clim_value);
%         xlim([-0.1 0.5])
%         xline(0,'LineStyle',':')
%         currentAx = gca; % 获取当前轴
%         subplotPosition = currentAx.Position; % 获取位置和大小
%         maxh=subplotPosition(4);
%         maxb=subplotPosition(2);
% 
%         subplotPosition(4)=maxh/max_num*size(used_plot_all_selected_1,1);
%         subplotPosition(2)=maxb+maxh-maxh/max_num*size(used_plot_all_selected_1,1);
%         ax.Position=subplotPosition;
%         axis off
% 
% 
%         if curr_stim==2
% 
%             cb = colorbar('south');  % 横向放在下方
%             cb.Position = [subplotPosition(1) subplotPosition(2)-0.05 0.2 0.01];
%             cb.Label.String = '\DeltaFR/FR_{0}';   % 给 colorbar 加标签
% 
%         end
%         if curr_stim==1 &curr_group==2
% 
%             pos = ax.Position;   % [x0 y0 width height]
% 
%             % 横线 (annotation 是 figure 坐标)
%             annotation('line', [pos(1)-0.05 pos(1)], [pos(2)-0.02 pos(2)-0.02], ...
%                 'Color','k','LineWidth',2);
%             % 竖线
%             annotation('line', [pos(1)-0.05 pos(1)-0.05], [pos(2)-0.02 pos(2)-0.02+100*maxh/max_num], ...
%                 'Color','k','LineWidth',2);
% 
%             % 文字
%             annotation('textbox', [pos(1)-0.1 pos(2)-0.07 0.2 0.05], ...
%                 'String','0.5 s', 'EdgeColor','none', ...
%                 'HorizontalAlignment','center');
%             annotation('textbox', [pos(1)-0.05 pos(2) 0.3 0.05], ...
%                 'String','100 neurons', 'EdgeColor','none', ...
%                 'HorizontalAlignment','center','Rotation',90);
% 
%         end
% 
% 
%         ax=subplot(4,3,6*curr_stim-6+3)
%         hold on
%         temp_plot=cell2mat(cellfun(@(x) nanmean(x,1)    ,used_filter_plot_all,'UniformOutput',false));
%         temp_mean=nanmean(temp_plot,1);
%         temp_error=std(temp_plot,0,1,'omitmissing')./sqrt(size(temp_plot,1));
%         ap.errorfill(t_bins,temp_mean,temp_error,colors{curr_group},0.5,0.1)
% 
%         ylim(yscale)
%         xlim([-0.1 0.5])
%         xline(0,'LineStyle',':')
% 
%         % axis off
% 
% 
% 
% 
%     end
% 
% 
%     means = [nanmean(proportion_response{1}, 1);nanmean(proportion_response{2}, 1)]';
%     sems = [std(proportion_response{1}, 0, 1,'omitmissing') ./ sqrt(size(proportion_response{1},1));...
%         std(proportion_response{2}, 0, 1,'omitmissing') ./ sqrt(size(proportion_response{2},1))]';
%     p =  ranksum(proportion_response{1}, proportion_response{2})
% 
%     ax=subplot(4,3,6*curr_stim)
%     ax.Color = 'none';    % 设置背景透明
% 
%     hold on
%     bar_handle = bar(1:2,means, 'grouped');
%     bar_handle.FaceColor = 'none';  % 允许每个柱子单独设色
%     bar_handle.EdgeColor = 'flat';  % 允许每个柱子单独设色
%     bar_handle.CData(1,:) = colors{1} ;  % 第一个柱子的颜色（RGB）
%     bar_handle.CData(2,:) = colors{2} ;  % 第二个柱子的颜色（RGB）
% 
%     errorbar(1:2, means, sems, 'k.', 'LineWidth', 1);
% 
%     % 添加散点
%     arrayfun(@(g) scatter(g*ones(length(proportion_response{g}),1) + randn(size(proportion_response{g},1),1)*0.05,...
%         proportion_response{g}, ...
%         20, 'filled', ...
%         'MarkerFaceColor', colors{g}), 1:2);
%     if p < 0.05
%         stars = repmat('*',1,sum(p<[0.05 0.01 0.001]));
%         y_sig = max(vertcat(proportion_response{:})) + 0.05;
%         plot(1:2, [1 1]*y_sig, 'k-');
%         text(1.5, y_sig+0.02, stars, 'HorizontalAlignment','center');
%     end
%     xticklabels({})
%     ylabel('fraction','FontWeight','normal')
%     ylim(bar_scale)
%     box off
%     drawnow
% end
% 
% exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s7a.eps'), ...
%     'ContentType','vector');
% 
% % clearvars('-except',main_preload_vars{:});
% 
% 
% %% group average by depth_version 1
% main_preload_vars = who;
% groups={'VA','AV'}
% titles={'L','M','V passive','4k','A passive','12k','V task','A task','iti move'};
% p_val=0.95
% for state=2
%     switch state
%         case 1
%             all_stim=[3 5];
%         case 2
%             all_stim=[7 8];
% 
%     end
% 
%     colors={[84 130 53]./255,[112  48 160]./255};
%     % colors={[0.3 0.3 1],[1 0.3 0.3]};
% 
%     colors1={[0.1706    0.1275    0.1165],[0.3294    0.5098    0.2078], [0.7451    0.8667    0.6706];...
%         [0.2706    0.0353    0.4667],[0.4392    0.1882    0.6275], [0.8196    0.7216    0.9019]};
%     colors_image={'G','P'}
%     z_min = 0;
%     z_max = 250;
%     bin_size_z = 10; % 单位：μm，根据实际尺度调整
% 
%     z_edges = [z_min:bin_size_z:z_max,inf];
% 
%     fig2 = figure('Position',[50 50 400 300]);
%     tl2 = tiledlayout(length(all_stim),3);
%     firing_rates_max_mice=cell(2,1);
%     for curr_group=1:2
% 
%         switch curr_group
%             case 1
%                 used_animals=animals(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
%                 used_animals_idx=anterior_learned_idx_VA(~cellfun(@isempty, anterior_learned_idx_VA','UniformOutput',true));
%             case 2
%                 used_animals=animals(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
%                 used_animals_idx=anterior_learned_idx_AV(~cellfun(@isempty, anterior_learned_idx_AV','UniformOutput',true));
%         end
% 
%         temp_single_idx=cell(length(used_animals),1);
%         temp_single_plot=cell(length(used_animals),1);
%         temp_probe_position=cell(length(used_animals),1);
%         temp_response=cell(length(used_animals),1);
%         temp_response_plot=cell(length(used_animals),1);
%         temp_cell_position=cell(length(used_animals),1);
%         for curr_animal=1:length(used_animals)
%             animal=used_animals{curr_animal};
%             temp_file_name=matfile([Path '\data\ephys_data\' animal '_ephys.mat']);
%             temp_single_plot{curr_animal}=temp_file_name.plot_single(used_animals_idx{curr_animal},1);
%             temp_single_idx{curr_animal}=temp_file_name.plot_idx(used_animals_idx{curr_animal},1);
%             temp_probe_position{curr_animal}=temp_file_name.all_celltypes(used_animals_idx{curr_animal},1);
%             temp_response{curr_animal}=temp_file_name.all_event_response_idx(used_animals_idx{curr_animal},1);
%             temp_response_plot{curr_animal}=temp_file_name.all_event_response_signle_neuron(used_animals_idx{curr_animal},1);
%             temp_cell_position{curr_animal}= cellfun(@(x,y) x-y , temp_file_name.all_cell_ccf_position_sorted(used_animals_idx{curr_animal},1),...
%                 temp_file_name.striatal_surface_position(used_animals_idx{curr_animal},1),'UniformOutput',false);
%         end
% 
%         single_neuron_each_rec_1=vertcat(temp_response_plot{:});
%         single_neuron_all_plot=cat(1,single_neuron_each_rec_1{:});
%         response_each_rec=vertcat(temp_response{:});
%         response_all=cat(1,response_each_rec{:});
%         single_neuron_each_position=vertcat(temp_cell_position{:});
%         single_neuron_position_all=cat(1,single_neuron_each_position{:});
% 
%         neuron_count_map_all=cell(2,1);
%         neuron_count_map=cell(2,1);
%         neuron_count_map_overlay=cell(2,1);
%         for curr_sorting=1:length(all_stim)
% 
% 
% 
%             used_stim=all_stim(curr_sorting);
% 
%             neuron_coords_all= cellfun(@(x)   x(:,2) ,single_neuron_each_position,'UniformOutput',false );
%             neuron_coords_each= cellfun(@(x,y) x(y(:,used_stim)>p_val,2),single_neuron_each_position,response_each_rec,'UniformOutput',false)
% 
% 
%             if used_stim==3|used_stim==5
%                 temp_single_idx=[3 5]
%             else
%                 temp_single_idx=[7 8]
%             end
% 
% 
% 
% 
%             firing_rates=cellfun(@(x) x(:,:,used_stim) ,single_neuron_each_rec_1,'UniformOutput',false);
% 
%             % === Step 1: 投影到冠状面（y-z） ===
%             projected_coords = neuron_coords_each; % 取 y 和 z
%             projected_coords_all = neuron_coords_all; % 取 y 和 z
% 
% 
%             % 分配每个神经元的 bin 索引
%             [neuron_count_map_all{curr_sorting},~,binIdx_all] = cellfun(@(x) histcounts(x, z_edges),projected_coords_all,'UniformOutput',false);
%             [neuron_count_map{curr_sorting},~,binIdx]= cellfun(@(x) histcounts(x, z_edges),projected_coords,'UniformOutput',false);
% 
% 
% 
%             firing_rates_bins =cellfun(@(x,idx) arrayfun(@(col) ...
%                 accumarray(idx, x(:,col), [length(z_edges)-1,1], @mean, NaN), ...
%                 1:size(x,2), 'UniformOutput', false),firing_rates,binIdx_all, 'UniformOutput', false);
% 
%             firing_rates_bins1=cellfun(@(x) cat(2,x{:}),firing_rates_bins,'UniformOutput',false);
%             firing_rates_bins2=nanmean(cat(3,firing_rates_bins1{:}),3);
% 
% 
%             figure(fig2)
%             ax2=nexttile(tl2,3*curr_sorting-3+curr_group)
%             % ax2=nexttile(tl2,curr_stim_now-3+3*curr_group)
% 
%             h2=imagesc(t_bins,z_edges(1:end-1),firing_rates_bins2)
%             xlim([-0.1 0.5]);
%             xticks([-0.1 0.5]);
%             hold on
%             xline(0,'LineStyle',':')
%             xlabel('time (s)');
%             ylim([z_edges(1)-0.5*bin_size_z  z_edges(end-1)+0.5*bin_size_z])
% 
%             if used_stim<7
%                 clim([0 2])
%             else
%                 clim([0 4])
%             end
%             colormap(ax2,ap.colormap(['W' colors_image{curr_group}]))
% 
%             % if curr_stim_now==2
%             %     colorbar('southoutside')
%             % end
% 
%             if curr_group==1
%                 % title(titles(used_stim),'FontWeight','normal')
%                 ylabel('depth (\mum)');
%                 yticks([z_edges(1) z_edges(end-1)]);
%             else
%                 yticks([]);
% 
%             end
% 
% 
%             % MUA max
%             ax1=nexttile(tl2,3*curr_sorting)
% 
%             firing_rates_max=cellfun(@(x) max(x(:,psth_use_t_stim),[],2),firing_rates_bins1,'UniformOutput',false );
%             firing_rates_max_mice{curr_sorting}{curr_group}=cat(2,firing_rates_max{:});
%             ap.errorfill(z_edges(1:end-1) , smoothdata(nanmean(cat(3,firing_rates_max{:}),3),'gaussian',4),...
%                 smoothdata( std(cat(3,firing_rates_max{:}),0,3,'omitmissing')./sqrt(size(cat(3,firing_rates_max{:}),3)),'gaussian',4),...
%                 colors{curr_group},0.1,0.5);
%             % ap.errorfill(z_edges(1:end-1) , nanmean(cat(2,firing_rates_max{:}),2),...
%             %    std(cat(3,firing_rates_max{:}),0,3,'omitmissing')./sqrt(size(cat(3,firing_rates_max{:}),3)),...
%             %   colors{curr_group},0.1,0.5);
% 
%             if used_stim<7
%                 ylim(ax1,[0 3])
%             else
%                 ylim(ax1,[0 8])
%             end
%             xlim(ax1,[z_edges(1) z_edges(end-1)])
% 
%             % title('activity','FontWeight','normal')
%             ylabel('\DeltaFR/FR_0')
%             xline(150,'LineStyle',':');
%             xline(50,'LineStyle',':')
%             xticks([]);
%             % set(gca, 'YAxisLocation', 'right')  % 把 x 轴移到上面
% 
%             view(ax1,90, 90);
% 
% 
%         end
% 
% 
%     end
% 
%     firing_rates_max_mice{1}
% 
%     cellfun(@(x) arrayfun(@(id) ds.shuffle_test(x{1}(id,:),x{2}(id,:),0,1),...
%         1:25,'UniformOutput',false),firing_rates_max_mice,'UniformOutput',false)
% 
%     cellfun(@(x) arrayfun(@(id) ds.shuffle_test(x{1}(id,:),x{2}(id,:),0,1),...
%         1:25,'UniformOutput',false),firing_rates_max_mice,'UniformOutput',false)
% 
% 
%     arrayfun(@(id) ds.shuffle_test(firing_rates_max_mice{1}{2}(id,:),firing_rates_max_mice{2}{2}(id,:),1,1),...
%         1:25,'UniformOutput',false)
%     %  switch state
%     %         case 1
%     % exportgraphics(gcf, fullfile(Path,'figures\eps\Fig 6d.eps'), ...
%     %     'ContentType','vector');
%     % case 2
%     % exportgraphics(gcf, fullfile(Path,'figures\eps\Fig s7b.eps'), ...
%     %     'ContentType','vector');
%     %     end
% end
% % clearvars('-except',main_preload_vars{:});

%% projection distribution
main_preload_vars = who;

allenAtlasPath = fileparts(which('template_volume_10um.npy'));
saveLocation=fullfile(allenAtlasPath,'temp_connect')
 % saveLocation = 'C:\Users\dsong\Documents\temp_connect'; % where to save the data downloaded from the Allen Connectivity dataset 
% allenAtlasPath =  'C:\Users\dsong\Documents\GitHub\osfstorage-archive'; % download from: https://figshare.com/articles/dataset/Modified_Allen_CCF_2017_for_cortex-lab_allenCCF/25365829 
fileName = ''; % leave empty to recompute each time (e.g. load the Allen raw data and sumnmarize it into one matrix), 
 
all_inputRegions={{'VIS'}, {'AUD'}}
corlors={'B','R'}
for curr_fig=1:2
% inputRegions = {'VIS'};
inputRegions = all_inputRegions{curr_fig};

mouseLine = ''; % leave empty to include all. use allen mouse line ids. 0 = wild-type. 
primaryInjection = true; % boolean, search for injections where 'injection' was the primary or not

experimentIDs = bsv.findConnectivityExperiments(inputRegions, mouseLine, primaryInjection);
% Fetch/load experiment data 
subtractOtherHemisphere = false;
loadAll = true; % if true, will load a 132 x 80 x 114 x number of experiments matrix instead of 132 x 80 x 114.
normalizationMethod = 'injectionVolume'; %  can be 'none', 'injectionIntensity' or 'injectionVolume'
groupingMethod = ' '; % leave empty or 'NaN' to average images all together. Other options include averaging by
% 'brainRegion', 'AP', 'ML', 'DV'

[experimentImgs, injectionSummary, experimentImgs_perExperiment] = bsv.fetchConnectivityData(experimentIDs, saveLocation, fileName, normalizationMethod,...
    subtractOtherHemisphere, groupingMethod, allenAtlasPath, loadAll);
% Plot projection data (in 2D) 
numberOfSlices =10; % for plotting purposes: divide target (output) structure into this many slices
numberOfPixels = 15; % for plotting purposes: divide each slice in target region in numberOfPixels x numberOfPixels
outputRegions = {'CP'}; % target region of interest


color=ap.colormap(['W' corlors{curr_fig}]);
plane = 'coronal'; % - not implemented yet - coronal or sagital
smoothing = 2; % - not implemented yet - none or a number (of pixels)
colorLimits = 'global'; % - not implemented yet - global, per slice or two numbers  
regionOnly = true; % - not implemented yet - whether to plot only one region or whole slices of the brain
% Plot!

bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color);

exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s5_' num2str(curr_fig)  '.eps']), ...
    'ContentType','vector');
end
 clearvars('-except',main_preload_vars{:});
