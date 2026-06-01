%%  EDF 5  different visual task types

main_preload_vars = who;

load(fullfile(Path,'data\revision\different_task_behavior.mat'));

colors={[0.4 0.4 1],[0.8 0.8 1],[0.5 0.5 1],[0.44 0.4 1],[0.44 0.4 1]};
use_group={[1 2],3}
for curr_groups=1:2
    figure('Position',[50 50 600 150]);
    tiledlayout(1,4)
    nexttile
    hold on
    for curr_group=use_group{curr_groups}
        cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.reaction_time{curr_group},'UniformOutput',false )

        cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.reaction_time{curr_group},...
            behavior_data.p_val{curr_group},'UniformOutput',false )
    end
    % ap.errorfill(1:length(RT_mean),RT_mean,RT_error,[0 0 0])
    xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
    yticks([0.1 0.5 1 2 4 6])
    set(gca, 'YScale', 'log')
    ylim([0.1 10])

    ylabel('Reaction time (s)')
    xlabel('Days')
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    nexttile
    hold on
    for curr_group=use_group{curr_groups}

        cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.performance{curr_group},'UniformOutput',false )

        cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.performance{curr_group},...
            behavior_data.p_val{curr_group},'UniformOutput',false )
    end

    xlim([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])
    ylim([-0.2 1])
    ylabel('Performance')
    xlabel('Days')
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])
    yticks([0 1])

    % sgtitle(title_name)


    temp_wf_task=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(3,a{end-1:end}),3)...
        ,x,'UniformOutput',false) ,behavior_data.task_kernels(use_group{curr_groups}),'UniformOutput',false)

    tem_image_task=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
        ,temp_wf_task,'UniformOutput',false);

    image_plot_task=feval(@(w)    nanmean(max(w(:,:,t_kernels>0&t_kernels<0.2 ,:),[],3),4) ,...
        feval(@(c) cat(4,c{:}), cat(2,tem_image_task{:})));

    nexttile
    imagesc(image_plot_task)
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));



    temp_wf_passive=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(4,a{end-1:end}),4)...
        ,x,'UniformOutput',false) ,behavior_data.passive_kernels(use_group{curr_groups}),'UniformOutput',false)

    tem_image_passive=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1),:),x),a,'UniformOutput',false)...
        ,temp_wf_passive,'UniformOutput',false);

    image_plot_passive=feval(@(w)    permute(nanmean(max(w(:,:,t_kernels>0&t_kernels<0.2 ,:,:),[],3),5),[1,2,4,3]) ,...
        feval(@(c) cat(5,c{:}), cat(2,tem_image_passive{:})));

    nexttile
    imagesc(image_plot_passive(:,:,3))
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));

    exportgraphics(gcf, fullfile(plab.locations.server_path,...
        ['Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_E' num2str(curr_groups) '.eps']), ...
        'ContentType','vector');
end

use_group=[ 4 5];
figure('Position',[50 50 1000 150]);
mainfig=tiledlayout(1,4)

plot_fig=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'none', 'Padding', 'none');
plot_fig.Layout.Tile = 1;  % 明确放在主 layout 的第 1 个 tile

for curr_group=use_group
    nexttile(plot_fig)
    hold on
    cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.reaction_time{curr_group},'UniformOutput',false )

    cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.reaction_time{curr_group},...
        behavior_data.p_val{curr_group},'UniformOutput',false )
    % ap.errorfill(1:length(RT_mean),RT_mean,RT_error,[0 0 0])
    set(gca, 'YScale', 'log')

    ylim([0.1 10])
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    switch curr_group
        case 4
            xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            yticks([0.1 0.5 1 5 10])
            ylabel('Reaction time (s)')
            xlabel('Days')
        case 5
            xlim([0 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            set(gca,'YColor','none')
    end
end

plot_fig2=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'none', 'Padding', 'none');
plot_fig2.Layout.Tile = 2;  % 明确放在主 layout 的第 1 个 tile
for curr_group=use_group
    nexttile (plot_fig2)
    hold on
    cellfun(@(x) plot(x,'Color',colors{curr_group},'LineWidth',2) ,behavior_data.performance{curr_group},'UniformOutput',false )
    cellfun(@(x,y) plot(find(y),x(y),'Color','r','LineWidth',2) ,behavior_data.performance{curr_group},...
        behavior_data.p_val{curr_group},'UniformOutput',false )

    ylim([-0.2 1])
    xticks([1 max(cellfun(@length ,behavior_data.performance{curr_group},'UniformOutput',true))])

    switch curr_group
        case 4
            xlim([1 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            yticks([0.1 0.5 1 2 4 6])
            % set(gca, 'YScale', 'log')
            ylabel('Performance')
            xlabel('Days')
        case 5
            xlim([0 max(cellfun(@length ,behavior_data.reaction_time{curr_group},'UniformOutput',true))])
            set(gca,'YColor','none')
    end
end
% sgtitle(title_name)


temp_wf_task=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(3,a{end-1:end}),3)...
    ,x,'UniformOutput',false) ,behavior_data.task_kernels(use_group),'UniformOutput',false)
tem_image_task=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
    ,temp_wf_task,'UniformOutput',false);
image_plot_task=cellfun(@(x) feval(@(v)  nanmean(v,3) , feval(@(c) cat(3,c{:}) ,cellfun(@(a)  max(a(:,:,t_kernels>0&t_kernels<0.2 ,:),[],3),x,'UniformOutput',false)))...
    ,tem_image_task,'UniformOutput',false  );

plot_fig3=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'tight', 'Padding', 'tight');
plot_fig3.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile

for curr_image=1:2
    nexttile (plot_fig3)
    imagesc(image_plot_task{curr_image})
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));
end


temp_wf_passive=   cellfun(@(x)   cellfun(@(a)      nanmean(cat(4,a{end-1:end}),4)...
    ,x,'UniformOutput',false) ,behavior_data.passive_kernels(use_group),'UniformOutput',false)
tem_image_passive=   cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1),:),x),a,'UniformOutput',false)...
    ,temp_wf_passive,'UniformOutput',false);
image_plot_passive=cellfun(@(x) feval(@(v)  nanmean(v,4) , feval(@(c) cat(4,c{:}) ,cellfun(@(a)  permute (max(a(:,:,t_kernels>0&t_kernels<0.2 ,:,:),[],3),[1,2,4,3]),x,'UniformOutput',false)))...
    ,tem_image_passive,'UniformOutput',false  );
plot_fig4=tiledlayout(mainfig,1, 2, ...
    'TileSpacing', 'tight', 'Padding', 'tight');
plot_fig4.Layout.Tile = 4;  % 明确放在主 layout 的第 1 个 tile

for curr_image=1:2
    nexttile (plot_fig4)
    imagesc(image_plot_passive{curr_image}(:,:,3))
    axis image off
    clim( 0.0003*[-1,1]);
    ap.wf_draw('ccf',[0.5 0.5 0.5]);
    colormap( ap.colormap(['KWB']));
end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_E3.eps'), ...
    'ContentType','vector');
 clearvars('-except',main_preload_vars{:});
