%% EDF 10 correlation between task vs passive kernels
main_preload_vars = who;
temp_data_all=cell(2,1);

for curr_group=1:2
    switch curr_group
        case 1
            animals = {'DS007','DS010','AP019','AP021','DS011','AP022'};
        case 2
            animals = {'DS000','DS004','DS014','DS015','DS016'};
    end

    temp_data_all{curr_group}=table;
    for curr_animal=1:length(animals)
        preload_vars=who;
        animal=animals{curr_animal};
        data_all=load(fullfile(Path,'data','revision','single_data',[animal '_all_data.mat']));

        switch curr_group
            case 1
                select_id_3=(strcmp([data_all.task_name],'stim_wheel_right_stage1')|strcmp([data_all.task_name],'stim_wheel_right_stage2'))&...
                    ~cellfun(@isempty ,data_all.wf_lcr_passive);
            case 2
                select_id_3=(strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume')|strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume'))&...
                    ~cellfun(@isempty ,data_all.wf_hml_passive_audio);


        end
        temp_p_val=arrayfun(@(id)  data_all.behavior_task{id}.rxn_l_p(1)<0.05, find(select_id_3),'UniformOutput',true);


        % wf_task

        all_wf_task=[data_all.wf_task(select_id_3)];
        % all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.iti_move_kernels{1}),all_wf_task,'UniformOutput',false);
        all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.stim_kernels{1}),all_wf_task,'UniformOutput',false);
        % all_task_image=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.stim_kernels{1},1)),x.all_iti_move_kernels{1}),all_wf_task,'UniformOutput',false);

        temp_data_all{curr_group}.task_image{curr_animal}=cat(4,all_task_image{find(temp_p_val==1,2,'last')});
        % wf_passive
        switch curr_group
            case 1
                all_wf_passive_1=[data_all.wf_lcr_passive(select_id_3)];
                used_passive=3;
            case 2
                all_wf_passive_1=[data_all.wf_hml_passive_audio(select_id_3)];
                used_passive=2;
        end
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x.kernels_decoding{used_passive},1)),x.kernels_decoding{used_passive}),all_wf_passive_1,'UniformOutput',false);

        temp_data_all{curr_group}.passive_image{curr_animal}=cat(4,all_passive_image_pre{find(temp_p_val==1,2,'last')})

        clearvars('-except',preload_vars{:});
    end

end

pairs={[1 1 ],[2 2],[1 2],[2 1]}
pairs_names={'Visual task vs Visual passive','Auditory task vs Auditory passive',...
    'Visual task vs Auditory passive','Auditory task vs Visual passive'};
temp_val=cell(4,1);
figure('Position',[50 50 600 900])
pair_colors={'B','R','K','K'}
for curr_pair=1:length(pairs)


    % task_images_max=feval(@(a)  cat(3,a{:}), ...
    %     cellfun(@(x)  permute(max(x(:,:,period,:),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(1)}.task_image,'UniformOutput',false));
    %
    task_images_max=feval(@(a)  cat(3,a{:}), ...
        cellfun(@(x)  permute(max(nanmean(x(:,:,kernels_period,:),4),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(1)}.task_image,'UniformOutput',false));


    task_length=size(task_images_max,3);

    % passive_images_max=feval(@(a)  cat(3,a{:}), ...
    %     cellfun(@(x)  permute(max(x(:,:,period,:),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(2)}.passive_image,'UniformOutput',false));
    %

    passive_images_max=feval(@(a)  cat(3,a{:}), ...
        cellfun(@(x)  permute(max(nanmean(x(:,:,kernels_period,:),4),[],3),[1,2,4,3])   ,temp_data_all{pairs{curr_pair}(2)}.passive_image,'UniformOutput',false));

    passive_length=size(passive_images_max,3);


    X=cat(3,task_images_max,passive_images_max);


    % ap.imscroll(X)
    % axis image off;
    % ap.wf_draw('ccf', [0.5 0.5 0.5]);
    % clim(0.0003 .* [ 0, 1]);
    % colormap( ap.colormap('WB' ));


    X = double(X);   % n*m*10

    K = size(X, 3);  % 这里应该是 10
    C = zeros(K, K);

    for i = 1:K
        vi = X(:,:,i);
        vi = vi(:);

        for j = 1:K
            vj = X(:,:,j);
            vj = vj(:);

            % 去掉常数向量导致的 NaN
            if std(vi) == 0 || std(vj) == 0
                C(i,j) = NaN;
            else
                C(i,j) = corr(vi, vj);
            end
        end
    end

    temp_val{curr_pair}=C(1:task_length,task_length+1:task_length+passive_length);


    K = size(C,1);

    % 创建mask：保留对角线 + 上三角
    mask = triu(true(K));

    % 其余设为 NaN（不显示）
    C_plot = nan(K);
    C_plot(mask) = C(mask);


   a1=nexttile
    imagesc(C_plot);
    colormap(a1, ap.colormap(['W' pair_colors{curr_pair}] ));
    set(gca, 'XAxisLocation', 'top');
    set(gca, 'YAxisLocation', 'right');
    axis image;
    box off
    clim([0 1])

    labels={'task','passive'}
    set(gca, 'XTick', [task_length/2 task_length+passive_length/2], 'YTick', [task_length/2 task_length+passive_length/2]);
    set(gca, 'XTickLabel', labels, 'YTickLabel', labels);
    % set(gca, 'TickLabelInterpreter', 'none');

    line([task_length+0.5 (task_length+passive_length)+0.5 ], [task_length+0.5 task_length+0.5 ],'Color',[1 0 0]);
    line([ task_length+0.5 task_length+0.5],[0.5 task_length+0.5 ],'Color',[1 0 0]);
    title(pairs_names{curr_pair})
end
colorbar;
caxis([0 1]);   % correlation 范围

temp_v=cellfun(@(x)  x(:),temp_val,'UniformOutput',false);
% figure('Position',[50 50 300 400])
nexttile
ds.make_bar_plot(cellfun(@(x)  x(:),temp_val,'UniformOutput',false),'ShowDots',0)
set(gca, 'XTick',[1:4],'XTickLabel', pairs_names);
ylim([0 1])
ylabel('Correlation')

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_C1.eps'), ...
    'ContentType','vector');

ds.shuffle_test(temp_v{3},temp_v{4})



% correlation trace
image_corr=cell(2,1);
for curr_group=1:2
    task_images=cat(4,temp_data_all{curr_group}.task_image{:});
    passive_images=cat(4,temp_data_all{curr_group}.passive_image{:});
    % image_corr=cell(size(task_images,4),1);
    for curr_image=1:size(task_images,4)
        used_task=reshape(task_images(:,:,:,curr_image),[],size(task_images,3))';
        used_passive=reshape(passive_images(:,:,:,curr_image),[],size(passive_images,3))';
        Cmat=zeros(size(used_task,2),1);
        for curr_pixel=1:size(used_task,2)
            Cmat(curr_pixel) = corr(used_task(:,curr_pixel),used_passive(:,curr_pixel));   % 注意转置！
        end
        image_corr{curr_group}{curr_image} = reshape(Cmat, size(task_images,1), size(task_images,2));
    end
end

figure('Position',[50 50 400 200]);
group_colors={'B','R'}
for curr_group=1:2
    a1=nexttile
    plot_image=nanmean(cat(3,image_corr{curr_group}{:}),3);
    % plot_image(plot_image<0.3)=0;
    imagesc(plot_image)
    axis image off;
    ap.wf_draw('ccf', [0.5 0.5 0.5]);
    clim([ 0, 1]);
    colormap(a1, ap.colormap(['W' group_colors{curr_group}]));
    colorbar

end

exportgraphics(gcf, fullfile(plab.locations.server_path,...
    'Lab\Papers\Song_2025\submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_EDF_C2.eps'), ...
    'ContentType','vector');

 clearvars('-except',main_preload_vars{:});
