%%     EDF 11  block test
main_preload_vars = who;

load(fullfile(Path,'data\wf_block_test.mat') )
animals =     { 'AP030','AP032','DS030','DS031','DS029'};
temp_image_max=cell(length(animals),1);
temp_plot_tace=cell(length(animals),1);
temp_image_trace=cell(length(animals),1);
for curr_animal=1:length(animals)

    temp_image{1}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task1(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{3}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task2(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{2}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive1(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{4}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive2(1),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));

    temp_image{5}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),...
        feval(@(c) c{1}.task_signle(end),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));
    temp_image{6}= feval(@(v) v{1} ,cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x(:,:,3)),...
        feval(@(c) c{1}.passive_signle(end),wf_px_kernels_all.wf_px_kernels(curr_animal)),'UniformOutput',false));

    temp_image_max{curr_animal}= feval(@(x)permute(max(x(:,:,t_kernels>0 & t_kernels<0.2,:),[],3),[1,2,4,3]),cat(4,temp_image{:}));
    temp_plot_tace{curr_animal}= ds.make_each_roi(cat(4,temp_image{:}), length(t_kernels),roi1);

    temp_image_trace{curr_animal}=temp_image;
end


% for curr_animal=[1 2 3 4 5]
%     ap.imscroll(cat(4,temp_image_trace{curr_animal}{:}),t_kernels)
%     axis image off;
%     ap.wf_draw('ccf', [0.5 0.5 0.5]);
%     clim(0.0003 .* [ 0, 1]);
%     colormap( ap.colormap('WG' ));
% end

figure('Position',[50 50 800 200]);
tiledlayout(1,5)
for curr_image=1:4
    nexttile
    imagesc( feval(@(c) c(:,:,curr_image),  nanmean(cat(4,temp_image_max{:}),4) ));
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim(0.0003 .* [ 0, 1]);
    colormap( ap.colormap('WG' ));


end
plot_trace_mean = feval(@(c) nanmean(c,4 ),cat(4,temp_plot_tace{:}));
plot_trace_error = feval(@(c) std(c,0,4,'omitmissing' )./ sqrt(size(c,4)),cat(4,temp_plot_tace{:}));

plot_max_mean=feval(@(c) permute(nanmean(max(c(:,t_kernels>0 & t_kernels<0.2,:,:),[],2),4),[1,3,2]),cat(4,temp_plot_tace{:}))
plot_max_error=feval(@(c) permute(std(max(c(:,t_kernels>0 & t_kernels<0.2,:,:),[],2),0,4,'omitmissing')./ sqrt(size(c,4)),[1,3,2]),cat(4,temp_plot_tace{:}))

colors={[1 0 0],[0 0 0],[1 0.5 0.5],[0.5 0.5 0.5]}

% hold on
nexttile
for curr_roi=[ 1 3 ]
    hold on
    ap.errorfill(1:4,plot_max_mean(curr_roi,1:4),plot_max_error(curr_roi,1:4))
end

ylim(0.0001*[0.8 2.2])
xlim([1 4])
xticklabels({'task 1','passive 1','task 2','passive 2'})
legend({'','mPFC','','aPFC'},'Box','off','Location','northeastoutside')


exportgraphics(gcf, fullfile(Path,...
    'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s11.eps'), ...
    'ContentType','vector');
 clearvars('-except',main_preload_vars{:});
