%% EDF 3   reaction time matched kernels
main_preload_vars = who;

% different brain regions roi
figure('Position',[50 50 200 200]);
merged_max=roi1(1).data.mask  +roi1(3).data.mask*2+...
    roi1(7).data.mask*3+roi1(9).data.mask*4+roi1(14).data.mask*5;
imagesc(merged_max)

axis image off;
cmap = [1 1 1
    0.40 0.60 0.85   % 蓝
    0.45 0.75 0.45   % 绿
    0.95 0.75 0.35   % 橙黄
    0.90 0.45 0.45   % 红
    0.65 0.50 0.85   % 紫
    ];

colormap(cmap);
h=ap.wf_draw('ccf', [0.5 0.5 0.5]);
% ap.wf_draw('grid')
% ap.wf_draw('area')
exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_A1.eps'), ...
    'ContentType','vector');


clearvars('-except',main_preload_vars{:});
%%

main_preload_vars = who;

load(fullfile(Path,'data','revision','wf_task_decoding_concatnate.mat'));

colors = [
    0.23 0.30 0.75
    0.40 0.50 0.85
    0.60 0.70 0.90
    0.80 0.85 0.95
    0.90 0.93 0.98
    ];

groups={'VA','AV'};
passive_workflows={'lcr_passive','hml_passive_audio'};
passive_id={3,2};
task_name={'visual position','audio volume'}
used_roi= [7 10];
cmap = [
    0.3 0.3 0.8;
    0.8 0.3 0.3;
    0.3 0.8 0.3;
    0.75 0.7 0.3;
    ];% 或者 lines(5), turbo(5)
% Colors = mat2cell(cmap, ones(size(temp_vel,1),1), 3);
Colors = mat2cell(cmap, ones(4,1), 3);
%
figure('Position',[50 50 800 600])
mainLayout = tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'none');
for curr_group=1:2

    temp_wf=wf_stim_kernels_concat.wf_kernels(ismember(wf_stim_kernels_concat.group,groups{curr_group}));
    tem_image= feval(@(c) cat(2,c{:}),  cellfun(@(a)   cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),a,'UniformOutput',false)...
        ,temp_wf,'UniformOutput',false));
    temp_stim2move=feval(@(c) cat(2,c{:}), wf_stim_kernels_concat.stim_to_move(ismember(wf_stim_kernels_concat.group,groups{curr_group})));
    temp_stim2move_mean= nanmean(cellfun(@mean,temp_stim2move,'UniformOutput',true),2);
    temp_stim2move_error= std(cellfun(@mean,temp_stim2move,'UniformOutput',true),0,2)./sqrt(size(temp_stim2move,2));
    temp_vel=feval(@(c) cat(2,c{:}), wf_stim_kernels_concat.wheel_velocity(ismember(wf_stim_kernels_concat.group,groups{curr_group})));
    temp_vel2= feval(@(c)  arrayfun(@(id)  cat(2,c{id,:}) ,1:size(temp_vel,1),'uni',false),cellfun(@(x)  median(x,1,'omitmissing')',temp_vel,'UniformOutput',false));
    temp_vel_mean=cellfun(@(x)  median(x,2,'omitmissing'),temp_vel2,'UniformOutput',false  );
    temp_vel_error=cellfun(@(x)  std(x,0,2)./sqrt(size(x,2)),temp_vel2,'UniformOutput',false  );

    tem_trace=cellfun(@(x)   ds.make_each_roi(x,t_kernels,roi1),tem_image,'UniformOutput',false);
    tem_trace_mean=arrayfun(@(id)    median(cat(3,tem_trace{id,:}),3,'omitmissing')  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_error=arrayfun(@(id)    std(cat(3,tem_trace{id,:}),0,3)./sqrt(size(tem_trace,2))  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_peak=arrayfun(@(id)   permute(max(cat(3,tem_trace{id,:}),[],2),[3,1,2])  ,1:size(temp_vel,1),'UniformOutput',false);
    tem_trace_peak_error=arrayfun(@(id)    std(max(cat(3,tem_trace{id,:}),[],2),0,3)./sqrt(size(tem_trace,2))  ,1:size(temp_vel,1),'UniformOutput',false);




    temp_wf_raw=feval(@(a)  cat(1,a{:}), cellfun(@(a) arrayfun(@(id)  nanmean(cat(3,a{id,:}),3) ,1:4,'uni',false)  ,...
        cellfun(@(x)  cat(2,x{:})         ,...
        wf_stim_kernels_concat.wf_raw(ismember(wf_stim_kernels_concat.group,groups{curr_group})),'uni',false),...
        'UniformOutput',false));

    tem_image_raw=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),temp_wf_raw,'UniformOutput',false);
    tem_trace_raw=cellfun(@(x)   ds.make_each_roi(x,t_kernels,roi1),tem_image_raw,'UniformOutput',false);
    tem_image_raw_mean= arrayfun(@(a) nanmean(cat(3,tem_trace_raw{:,a}),3) ,1:4,'uni',false);
    tem_image_raw_error=arrayfun(@(a) nanstd(cat(3,tem_trace_raw{:,a}),0,3)./sqrt(size(tem_trace_raw,1)) ,1:4,'uni',false);
    tem_image_nanmean=feval(@(a)   cat(4,a{:}), arrayfun(@(a) nanmean(cat(4,tem_image_raw{:,a}),4) ,1:4,'uni',false));



    % tem_image_raw_mean= feval(@(b)   cat(4,b{:}), arrayfun(@(a) nanmean(cat(4,tem_image_raw{:,a}),4) ,1:4,'uni',false));
    % ap.imscroll(tem_image_raw_mean,t_task)
    % axis image off
    % clim( 0.02*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['KWB']));

    % figure('Position',[50 50 200 600]);
    imageLayout= tiledlayout(mainLayout,5,1,'Padding','tight','TileSpacing','none','TileIndexing','rowmajor')
    imageLayout.Layout.Tile = 2*curr_group-1;  % 明确放在主 layout 的第 1 个 tile

    for curr_roi=[used_roi(curr_group) 1 3  14]
        nexttile(imageLayout)
        for curr_state=2:size(temp_vel,1)
            hold on
            ap.errorfill(t_task,tem_image_raw_mean{curr_state}(curr_roi,:),...
                tem_image_raw_error{curr_state}(curr_roi,:),Colors{curr_state},0.5);
            xlim([-0.1 1])
            % if curr_roi

            if curr_roi==7 ||curr_roi==10
                temp_max_x=feval(@(a) max(a(curr_roi,:,:),[],'all'),cat(5,tem_image_raw_mean{:}))*1.2;
                temp_max = ceil(temp_max_x/10^floor(log10(abs(temp_max_x)))) * 10^floor(log10(abs(temp_max_x)));

                temp_min_x=feval(@(a) min(a(curr_roi,:,:),[],'all'),cat(5,tem_image_raw_mean{:}))*1.5;
                temp_min = ceil(temp_min_x/10^floor(log10(abs(temp_min_x)))) * 10^floor(log10(abs(temp_min_x)));
            else
                temp_max=0.015;
                temp_min=-0.005;
            end

            ylim([temp_min temp_max ])
            % ylim([-0.0002 0.0003])
            xline(0,'.k')
            xline(temp_stim2move_mean(curr_state),'Color',Colors{curr_state})
            axis off
        end
        line ([-0.1 0],[temp_min temp_min],'Color','k')
        line ([-0.1 -0.1],[temp_min 0],'Color','k')
        text(-0.15, temp_min*1.2,num2str(abs(temp_min)),'Rotation',90)
        text(-0.1, temp_min*1.2, '0.1 s')
    end

    nexttile(imageLayout)
    hold on
    cellfun(@(x,y,z) ap.errorfill(surround_time_points,x,y,z,0.5),...
        temp_vel_mean(2:4),temp_vel_error(2:4),Colors(2:4)','UniformOutput',false)
    temp_max2 = feval(@(a)  max(a,[],'all') ,  cat(1,temp_vel_mean{:}));
    temp_min2 = feval(@(a)  min(a,[],'all') ,  cat(1,temp_vel_mean{:}))*1.2;

    ylim([temp_min2 temp_max2])
    % ylim([-1500 2000])
    xlim([-0.1 1])
    xline(temp_stim2move_mean(curr_state),'.r')
    xline(0,'.k')
    axis off
    line ([-0.1 0],[-1000 -1000],'Color','k')
    line ([-0.1 -0.1],[-1000 -500],'Color','k')
    text(-0.15, -1000*1.2,num2str(abs(-500)),'Rotation',90)
    text(-0.1, -1000*1.2, '0.1 s')

    imageLayout=tiledlayout(mainLayout,5,1,'Padding','tight','TileSpacing','none','TileIndexing','rowmajor')

    imageLayout.Layout.Tile = 2*curr_group;  % 明确放在主 layout 的第 1 个 tile

    for curr_roi=[used_roi(curr_group) 1 3  14]
        nexttile(imageLayout)
        for curr_state=2:size(temp_vel,1)
            hold on
            ap.errorfill(t_kernels,tem_trace_mean{curr_state}(curr_roi,:),...
                tem_trace_error{curr_state}(curr_roi,:),Colors{curr_state},0.1);
            xlim([-0.1 1])
            % if curr_roi

            if curr_roi==7 ||curr_roi==10
                temp_max_x=feval(@(a) max(a(curr_roi,:,:),[],'all'),cat(5,tem_trace_mean{:}))*1.2;
                temp_min_x=feval(@(a) min(a(curr_roi,:,:),[],'all'),cat(5,tem_trace_mean{:}))*1.5;

                temp_max = ceil(temp_max_x/10^floor(log10(abs(temp_max_x)))) * 10^floor(log10(abs(temp_max_x)));

                temp_min = ceil(temp_min_x/10^floor(log10(abs(temp_min_x)))) * 10^floor(log10(abs(temp_min_x)));

            else
                temp_max=0.0004;
                temp_min=-0.0001;
            end

            ylim([temp_min temp_max ])
            % ylim([-0.0002 0.0003])
            xline(0,'.k')
            xline(temp_stim2move_mean(curr_state),'Color',Colors{curr_state})
            axis off
        end



        if curr_group==2 &curr_roi==used_roi(curr_group)
            h = findobj(gca,'Type','line');
            legend(h([1 2 3]), {'>0.3s','0.2-0.3s','0-0.2s'},'Box','off')
        end

        line ([-0.1 0],[temp_min temp_min],'Color','k')
        line ([-0.1 -0.1],[temp_min 0],'Color','k')
        text(-0.15, temp_min*1.2,num2str(abs(temp_min)),'Rotation',90)
        text(-0.1, temp_min*1.2, '0.1 s')

    end

    nexttile(imageLayout)
    hold on
    cellfun(@(x,y,z) ap.errorfill(surround_time_points,x,y,z,0.5),...
        temp_vel_mean(2:4),temp_vel_error(2:4),Colors(2:4)','UniformOutput',false)
    temp_max2 = feval(@(a)  max(a,[],'all') ,  cat(1,temp_vel_mean{:}));
    temp_min2 = feval(@(a)  min(a,[],'all') ,  cat(1,temp_vel_mean{:}))*1.2;


    ylim([temp_min2 temp_max2])
    % ylim([-1500 2000])
    xlim([-0.1 1])
    xline(temp_stim2move_mean(curr_state),'.r')
    xline(0,'.k')
    axis off
    line ([-0.1 0],[-1000 -1000],'Color','k')
    line ([-0.1 -0.1],[-1000 -500],'Color','k')
    text(-0.15, -1000*1.2,num2str(abs(-500)),'Rotation',90)
    text(-0.1, -1000*1.2, '0.1 s')



end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_A2.eps'), ...
    'ContentType','vector');

clearvars('-except',main_preload_vars{:});
