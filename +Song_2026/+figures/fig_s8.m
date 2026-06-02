%% EDF 8  visual passive size 20 60

 main_preload_vars = who;

load(fullfile(Path,'data','visual_size_passive_compare.mat'));

tem_passive_s_image=cell(2,1);
tem_passive_l_image=cell(2,1);
trace_s_mean=cell(2,1);
trace_l_mean=cell(2,1);
for curr_animal=1:5

    tem_passive_s=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),passive_data.lcr_passive(curr_animal),'UniformOutput',false);
    tem_passive_l=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),passive_data.lcr_passive_size60(curr_animal),'UniformOutput',false);
    tem_passive_s_image{curr_animal}=permute(nanmean(max(tem_passive_s{1}(:,:,kernels_period,3,:),[],3),5),[1,2,5,3,4]);
    tem_passive_l_image{curr_animal}=permute(nanmean(max(tem_passive_l{1}(:,:,kernels_period,3,:),[],3),5),[1,2,5,3,4]);

    trace_s=ds.make_each_roi(permute(tem_passive_s{1}(:,:,:,3,:),[1,2,3,5,4]),t_kernels,roi1);
    trace_l=ds.make_each_roi(permute(tem_passive_l{1}(:,:,:,3,:),[1,2,3,5,4]),t_kernels,roi1);
    trace_s_mean{curr_animal}=nanmean(trace_s,3);
    trace_s_error=std(trace_s,0,3,'omitmissing')./sqrt(size(trace_s,3));
    trace_l_mean{curr_animal}=nanmean(trace_l,3);
    trace_l_error=std(trace_l,0,3,'omitmissing')./sqrt(size(trace_l,3));

    % 
    % figure('Position',[50 50 800 200]);
    % tiledlayout(1,4)
    % 
    % nexttile
    % imagesc(tem_passive_s_image{curr_animal})
    % axis image off
    % clim( 0.0003*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['PWG']));
    % nexttile
    % imagesc(tem_passive_l_image{curr_animal})
    % axis image off
    % clim( 0.0003*[-1,1]);
    % ap.wf_draw('ccf',[0.5 0.5 0.5]);
    % colormap( ap.colormap(['PWG']));
    % 
    % nexttile
    % hold on
    % ap.errorfill(t_kernels,trace_s_mean{curr_animal}(7,:),trace_s_error(1,:))
    % ap.errorfill(t_kernels,trace_l_mean{curr_animal}(7,:),trace_l_error(1,:))
    % title('V1')
    % xlim([-0.1 0.5])
    % 
    % nexttile
    % hold on
    % ap.errorfill(t_kernels,trace_s_mean{curr_animal}(1,:),trace_s_error(1,:))
    % ap.errorfill(t_kernels,trace_l_mean{curr_animal}(1,:),trace_l_error(1,:))
    % title('mPFC')
    % xlim([-0.1 0.5])


end

trace_s_all_mean=nanmean(cat(3,trace_s_mean{:}),3);
trace_s_all_error=std(cat(3,trace_s_mean{:}),0,3)./sqrt(length(trace_s_mean));
trace_l_all_mean=nanmean(cat(3,trace_l_mean{:}),3);
trace_l_all_error=std(cat(3,trace_l_mean{:}),0,3)./sqrt(length(trace_l_mean));

tem_passive_l_image_all=nanmean(cat(3,tem_passive_l_image{:}),3);
tem_passive_s_image_all=nanmean(cat(3,tem_passive_s_image{:}),3);

max_data_s=feval(@(a) cat(2,a{:}), cellfun(@(x) max(x([1 7],kernels_period),[],2),   trace_s_mean,'UniformOutput',false))
max_data_l=feval(@(a) cat(2,a{:}), cellfun(@(x) max(x([1 7],kernels_period),[],2),   trace_l_mean,'UniformOutput',false))

% figure
% hold on
% scatter(max_data_s(2,:), max_data_s(1,:), 50,'filled')
% scatter(max_data_l(2,:), max_data_l(1,:), 50,'filled')
% xlim([0 0.001])
% ylim([0 0.001])

figure('Position',[50 50 800 200]);
tiledlayout(1,4)

nexttile
imagesc(tem_passive_s_image_all)
axis image off
clim( 0.0003*[-1,1]);
ap.wf_draw('ccf',[0.5 0.5 0.5]);
title('size 20')

colormap( ap.colormap(['PWG']));
nexttile
imagesc(tem_passive_l_image_all)
axis image off
clim( 0.0003*[-1,1]);
ap.wf_draw('ccf',[0.5 0.5 0.5]);
colormap( ap.colormap(['PWG']));
title('size 60')

nexttile
hold on
ap.errorfill(t_kernels,trace_s_all_mean(7,:),trace_s_all_error(7,:))
ap.errorfill(t_kernels,trace_l_all_mean(7,:),trace_l_all_error(7,:))
title('V1')
xlim([-0.1 0.5])
set(gca,'Color','none')

nexttile
hold on
ap.errorfill(t_kernels,trace_s_all_mean(1,:),trace_s_all_error(1,:))
ap.errorfill(t_kernels,trace_l_all_mean(1,:),trace_l_all_error(1,:))
title('mPFC')
xlim([-0.1 0.5])
legend({'','Small','','Large'},'Box','off','Location','northeastoutside')
set(gca,'Color','none')


exportgraphics(gcf, fullfile(Path,...
    'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s8.eps'))
 clearvars('-except',main_preload_vars{:});
