%%    EDF 16    mPFC aPFC ephys 
main_preload_vars = who;

animals={'DS029','DS030','DS031'};
temp_kernels=cell(length(animals),1);

for curr_animal=1:length(animals)
    animal=animals{curr_animal};
    temp_path=matfile(fullfile(Path,'data','single_data',[animal '_all_data.mat']));
    temp_name_idx=cellfun(@(x) strcmp(x,'stim_wheel_right_stage2_mixed_VA_earphone'),temp_path.task_name,'uni',true );
    temp_task=temp_path.wf_task;
    temp_kernels{curr_animal}= temp_task{~cellfun(@isempty ,temp_task,'UniformOutput',true)&temp_name_idx}.stim_kernels;
end

wf_task=...
    feval(@(aa,bb)    cat(4,aa,bb),...
    feval(@(c)  cat(3,c{:}) ,arrayfun(@(id)  temp_kernels{id}{1}   , 1:length(animals),'UniformOutput',false)),...
    feval(@(c)  cat(3,c{:}) ,arrayfun(@(id)  temp_kernels{id}{2}   , 1:length(animals),'UniformOutput',false)));

temp_image=plab.wf.svd2px(U_master(:,:,1:size(wf_task,1)),wf_task);
image_max=permute(nanmean(max(temp_image(:,:,kernels_period,:,:),[],3),4),[1,2,5,4,3]);

figure('Position',[50 50 300 600])
Color={'B','R'}
probe_colors={[0.5 1 0.5],[0.2 0.5 0.2],[0 0.8 0]}
h=struct
for curr_mod=1:2
    ax=nexttile
    imagesc(ax,image_max(:,:,curr_mod))
    axis image off;
    clim(0.0003 .* [0, 1]);
    colormap(ax, ap.colormap(['W' Color{curr_mod}] ));
    hold on
    for curr_animal=1:length(animals)

        ap.wf_draw('probe_insertion',animals{curr_animal});
        h.(['l' num2str(curr_animal)]) = findobj(gca,'Type','line');
        for i = 1:length( h.(['l' num2str(curr_animal)]))
            if strcmp(h.(['l' num2str(curr_animal)])(i).Marker,'.')&h.(['l' num2str(curr_animal)])(i).MarkerSize==20
                h.(['l' num2str(curr_animal)])(i).Color = probe_colors{curr_animal};
                h.(['l' num2str(curr_animal)])(i).MarkerSize = 10;
            end
        end
    end
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
end

exportgraphics(gcf, fullfile(Path,...
    ['submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s16a0.eps']), ...
    'ContentType','vector');

raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');
period=t_bins>0&t_bins<0.2;

animals={'DS029','DS030','DS031'};

% mPFC_idx={[2,4],[2,4,5],[3,4,5]}
responsive_visual_all=cell(2,1);
responsive_audio_all=cell(2,1);
responsive_or_all=cell(2,1);
visual_passive_all=cell(2,1);
audio_passive_all=cell(2,1);
visual_passive_all_2=cell(2,1);
audio_passive_all_2=cell(2,1);
for curr_group=1:2
    switch curr_group
        case 1
            PFC_idx={[1 ],[3],[1 2]};
        case 2
            PFC_idx={[2,4],[2,4,5],[3,4,5]};
    end


    % visual_task=cell(length(animals),1);
    % audio_task=cell(length(animals),1);
    visual_passive_1=cell(length(animals),1);
    audio_passive_1=cell(length(animals),1);
    visual_passive_2=cell(length(animals),1);
    audio_passive_2=cell(length(animals),1);

    responsive_idx=cell(length(animals),1);
    responsive_visual=cell(length(animals),1);
    responsive_audio=cell(length(animals),1);
    temp_ephys=cell(length(animals),1);
    for curr_animal=1:length(animals)
        animal=animals{curr_animal};
        temp_path=matfile(fullfile(Path,'data','single_data',[animal '_all_data.mat']));
        % temp_path.data_all_index


        temp_audio_passive=temp_path.ephys_hml_passive_audio_earphone;
        temp_visual_passive=temp_path.ephys_lcr_passive;
        temp_task=temp_path.ephys_task;
        temp_behavior=temp_path.behavior_task;
        ephys_idx=~cellfun(@isempty ,temp_visual_passive );

        ephys_audio=temp_audio_passive(ephys_idx);
        ephys_visual=temp_visual_passive(ephys_idx);
        ephys_task=temp_task(ephys_idx);
        behavior=temp_behavior(ephys_idx);




        for curr_day=1:length(PFC_idx{curr_animal})
            temp_day=PFC_idx{curr_animal}(curr_day);


            response_idx=ephys_visual{temp_day}.response_p{3}>0.95|ephys_audio{temp_day}.response_p{2}>0.95;
            depth_idx=ephys_task{temp_day}.depth>1840;
            responsive_idx{curr_animal}{curr_day}=response_idx(find(depth_idx));
            responsive_visual{curr_animal}{curr_day}=feval(@(a) a((find(depth_idx))), ephys_visual{temp_day}.response_p{3}>0.95);
            responsive_audio{curr_animal}{curr_day}=feval(@(a) a((find(depth_idx))), ephys_audio{temp_day}.response_p{2}>0.95);

            % visual_task{curr_animal}{curr_day}= ephys_task{temp_day}.psth(depth_idx,:,1);
            % audio_task{curr_animal}{curr_day}= ephys_task{temp_day}.psth(depth_idx,:,2);
            visual_passive_1{curr_animal}{curr_day}=ephys_visual{temp_day}.psth(depth_idx,:,3);
            audio_passive_1{curr_animal}{curr_day}=ephys_audio{temp_day}.psth(depth_idx,:,2);

            visual_passive_2{curr_animal}{curr_day}=ephys_visual{temp_day}.psth(depth_idx,:,9);
            audio_passive_2{curr_animal}{curr_day}=ephys_audio{temp_day}.psth(depth_idx,:,8);


            temp_ephys{curr_animal}{curr_day}= cat(3,...
                ephys_visual{temp_day}.psth(depth_idx,:,:),...
                ephys_audio{temp_day}.psth(depth_idx,:,:));

        end


    end


    responsive_visual_all{curr_group}=cat(2,responsive_visual{:});
    responsive_audio_all{curr_group}=cat(2,responsive_audio{:});
    responsive_or_all{curr_group}=cat(2,responsive_idx{:});
    visual_passive_all{curr_group}=cat(2,visual_passive_1{:});
    audio_passive_all{curr_group}=cat(2,audio_passive_1{:});

    visual_passive_all_2{curr_group}=cat(2,visual_passive_1{:});
    audio_passive_all_2{curr_group}=cat(2,audio_passive_1{:});
    % visual_task_all{curr_group}=cat(2,visual_task{:});
    % audio_task_all{curr_group}=cat(2,audio_task{:});
    % ephys_data_all{curr_group}=cat(2,temp_ephys{:});

end

temp_or_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_or_all,'uni',false);
temp_v_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_visual_all,'uni',false);
temp_a_idx=cellfun(@(x)  cat(1,x{:})  ,responsive_audio_all,'uni',false);
temp_v_frac=cellfun(@(x) cellfun(@(a) sum(a)/length(a) ,x,'UniformOutput',true) ,responsive_visual_all,'UniformOutput',false   );
temp_a_frac=cellfun(@(x) cellfun(@(a) sum(a)/length(a) ,x,'UniformOutput',true) ,responsive_audio_all,'UniformOutput',false   );
temp_frac=cellfun(@(x,y)    {x,y}, temp_v_frac,temp_a_frac,'UniformOutput',false   );

idx_A_only = cellfun(@(x,y)  find(x&~y),temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_AB     = cellfun(@(x,y) find (x&y),temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_B_only = cellfun(@(x,y)  find(~x&y),temp_v_idx,temp_a_idx,'UniformOutput',false);

temp_v_p=cellfun(@(x)     cat(1,x{:}), visual_passive_all,'UniformOutput',false);
temp_a_p=cellfun(@(x)     cat(1,x{:}), audio_passive_all,'UniformOutput',false);
temp_v_p_2=cellfun(@(x)     cat(1,x{:}), visual_passive_all_2,'UniformOutput',false);
temp_a_p_2=cellfun(@(x)     cat(1,x{:}), audio_passive_all_2,'UniformOutput',false);



[~,max_idx_v]=cellfun(@(x)    max(x(:,period) ,[],2)  ,temp_v_p_2,'UniformOutput',false);
[~,sort_idx_v] = cellfun(@(x) sortrows( x,"ascend"),max_idx_v,'uni',false);
[~,max_idx_a]=cellfun(@(x)    max(x(:,period) ,[],2)  ,temp_a_p_2,'UniformOutput',false);
[~,sort_idx_a] = cellfun(@(x) sortrows( x,"ascend"),max_idx_a,'uni',false);

idx_order_sort=cellfun(@(x,y,z,a) [ z(ismember(z,(find(x&~y))));...
    z(ismember(z,(find(x&y))));...
    a(ismember(a,(find(~x&y))))],...
    temp_v_idx,temp_a_idx,sort_idx_v,sort_idx_a,'UniformOutput',false);


 idx_order=cellfun(@(x,y)  [find(x&~y); find(x&y); find(~x&y)],temp_v_idx,temp_a_idx,'UniformOutput',false);
idx_orders=cellfun(@(x,y)  {find(x&~y); find(x&y); find(~x&y)},temp_v_idx,temp_a_idx,'UniformOutput',false);


temp_v_p_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    visual_passive_all,responsive_or_all,'UniformOutput',false);
temp_a_p_mean=cellfun(@(x,y) feval(@(m) cat(1,m{:}),  cellfun(@(a,b)  nanmean(a(b,:),1),x,y,'UniformOutput',false)) ,...
    audio_passive_all,responsive_or_all,'UniformOutput',false);


for curr_group=1:2
    figure;
    mainfig=tiledlayout(1,3,'TileSpacing','tight')

    A1=nexttile
    imagesc(t_bins,[],temp_v_p{curr_group}(idx_order_sort{curr_group},:) )
    % imagesc(t_bins,[],temp_v_p{curr_group}(idx_order{curr_group},:) )
    clim([0 3])
    xlim([-0.1 0.5])
    colormap(A1,ap.colormap(['WB']))
    yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    xline(0)
    colorbar('southoutside')
    axis off

    A2=nexttile
    imagesc(t_bins,[],temp_a_p{curr_group}(idx_order_sort{curr_group},:) )
        % imagesc(t_bins,[],temp_a_p{curr_group}(idx_order{curr_group},:) )

    clim([0 3])
    xlim([-0.1 0.5])
    colormap(A2,ap.colormap(['WR']))
    yline(size(idx_A_only{curr_group},1),'LineWidth',1,'Color',[0 0 0])
    yline(size([idx_A_only{curr_group}; idx_AB{curr_group}],1),'LineWidth',1,'Color',[0 0 0])
    xline(0)
    line ([-0.1 0],[1 1],'linewidth',1,'color','k')
    line ([-0.1 -0.1],[0 100],'linewidth',1,'color','k')
    text(-0.15, 150,'100 neurons','Rotation',90)
    text(-0.1, -10, '0.1 s')

    colorbar('southoutside')
    axis off


    plot_fig=tiledlayout(mainfig,2 ,1, ...
        'TileSpacing', 'tight');
    plot_fig.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile
    a4=nexttile(plot_fig,1)
    hold on
    ap.errorfill(t_bins,nanmean(temp_v_p_mean{curr_group},1),nanstd(temp_v_p_mean{curr_group},0,1)./sqrt(size(temp_v_p_mean{curr_group},1)),[0 0 1])
    ap.errorfill(t_bins,nanmean(temp_a_p_mean{curr_group},1),nanstd(temp_a_p_mean{curr_group},0,1)./sqrt(size(temp_a_p_mean{curr_group},1)),[1 0 0])
    ylim([-0.1 2.5])
    xlim([-0.1 0.5])
    xline(0)
    line([-0.1 -0.1],[0 0.5],'linewidth',1,'color','k')
    line([-0.1 0],[0 0],'linewidth',1,'color','k')
    text(-0.15, 0,'0.5','Rotation',90)
    text(-0.1, -0.2, '0.1 s')

    axis off

    a4=nexttile(plot_fig,2)
    ds.make_bar_plot(temp_frac{curr_group},'ColorCell',{[0 0 1],[1 0 0]})
    ylabel('Fraction')
    ylim([0 0.3])
    yticks([0 0.3])
    xticks([])
    set(gca,'Color','none')

    ds.shuffle_test(temp_frac{curr_group}{1},temp_frac{curr_group}{2},1)

    exportgraphics(gcf, fullfile(plab.locations.server_path,...
        ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_Ja' num2str(curr_group)   '.eps']), ...
        'ContentType','vector');

end

response_sort_both=...
    cellfun(@(a,b,c)     cellfun(@(x,y,z)   [x(z) y(z)] ,a,b,c,'uni',false)',...
    responsive_visual_all,responsive_audio_all,...
    responsive_or_all,'UniformOutput',false);

response_sort_all_both=cellfun(@(x)   cat(1,x{:}),response_sort_both,'UniformOutput',false);
response_both_fration_both = cellfun(@(a)  cellfun(@(x) sum(sum(x,2)==2)./size(x,1),a,'UniformOutput',true) ,...
    response_sort_both,'UniformOutput',false);

n_shuff=1000;
response_both_shuff_both= cellfun(@(x)  arrayfun(@(shuff) ap.shake(x,1),1:n_shuff,...
    'UniformOutput',false),response_sort_all_both,'UniformOutput',false);


response_both_fration_shuff_both=cellfun(@(a)  cellfun(@(x) sum( sum(x,2)==2)/ length(x),a,'UniformOutput',true),...
    response_both_shuff_both,'UniformOutput',false);
line_mean_both=cellfun(@(a)  (prctile(a,95)+ ...
    prctile(a,5))/2,response_both_fration_shuff_both,'UniformOutput',false);
line_error_both=cellfun(@(a) (prctile(a,95)-...
    prctile(a,5))/2,response_both_fration_shuff_both,'UniformOutput',false);







colors = {[0 0 1],[ 1.0, 0.647, 0.0],[1 0 0]};
title_names={'mPFC','aPFC'}

scales_all=[0 20 80];
scales = {[0 20 80],[0 8 20],[0 4 10]};
colors = {[0 0 1],[ 1.0, 0.647, 0.0],[1 0 0]};
for curr_group=1:2
    figure('Position',[50 50 200 200])
    ax1 = axes('Position',[0.2 0.2 0.5 0.5]); % 左边大一些
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    xlim(scales_all(1:2))
    ylim(scales_all(1:2))
    xticks(scales_all(1:2))
    yticks(scales_all(1:2))
    set(gca,'Color','none')
    xlabel({'visual response' ;'(\DeltaFR/FR_{0})'})
    ylabel({'auditory response';' (\DeltaFR/FR_{0})'})

    ax2 = axes('Position',[0.75 0.2 0.2 0.5]);
    hold on

    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    xlim(scales_all(2:3))
    ylim(scales_all(1:2))
    xticks(scales_all(3))
    set(gca,'YTick',[])
    set(gca,'YColor','none')
    set(gca,'Color','none')

    ax3 = axes('Position',[0.2 0.75 0.5 0.2]);
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    ylim(scales_all(2:3))
    xlim(scales_all(1:2))
    set(gca,'Color','none')
    yticks(scales_all(3))
    set(gca,'XTick',[])
    set(gca,'XColor','none')

    ax4 = axes('Position',[0.75 0.75 0.2 0.2]);
    hold on
    arrayfun(@(id) scatter(max(temp_v_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        max(temp_a_p{curr_group}(idx_orders{curr_group}{id},period),[],2),...
        20,'filled','MarkerFaceColor',colors{id},'MarkerFaceAlpha',0.5),[ 3 2 1],'uni',false)

    ylim(scales_all(2:3))
    xlim(scales_all(2:3))
    axis off
    title(title_names{curr_group},'FontWeight','normal')

    annotation('line',[0.74 0.76],[0.19 0.21],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.69 0.71],[0.19 0.21],'Color','k','LineWidth',1) % 右斜杠
    annotation('line',[0.19 0.21],[0.74 0.76],'Color','k','LineWidth',1) % 左斜杠
    annotation('line',[0.19 0.21],[0.69 0.71],'Color','k','LineWidth',1) % 右斜杠

exportgraphics(gcf, fullfile(Path,...
    ['submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s16f' num2str(curr_group)   '.eps']), ...
    'ContentType','vector'); 

end



figure('Position',[50 50 100 200])
colors={[0 0.5 0],[0.5 0 0.2]};
maker_size=20;

hold on
arrayfun(@(curr_group)  ap.errorfill(curr_group-0.4:0.8:curr_group+0.4,[line_mean_both{curr_group} line_mean_both{curr_group}],...
    [line_error_both{curr_group} line_error_both{curr_group}],[0.5 0.5 0.5],1),1:2);
ds.make_bar_plot(response_both_fration_both,'ColorCell',colors,'BarAlpha',0.2,'DotSize',maker_size)

ylim([0 0.5])
yticks([0 0.5])
xticks([1 2])
xticklabels({'mPFC','aPFC'})
ylabel(['(V∩A)/(V∪A)' ])
set(gca,'color','none')

exportgraphics(gcf, fullfile(Path,...
    ['submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s16g' num2str(curr_group)   '.eps']), ...
    'ContentType','vector'); 

 clearvars('-except',main_preload_vars{:});
