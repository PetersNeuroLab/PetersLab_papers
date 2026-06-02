%% EDF 14 nose movement in passive
main_preload_vars = who;

groups_name={'VA','AV'};
modes_name={'Visual','Auditory'};

all_data=struct;

for curr_group=1:2
    switch curr_group
        case 1
            animals = {'DS007','DS010','AP019','AP021','DS011','AP022'};
        case 2
            animals = {'DS000','DS004','DS015','DS016'};
    end
    for curr_mode=1:2

        temp_data_all=table;
        for curr_animal=1:length(animals)
            preload_vars=who;
            animal=animals{curr_animal};
            data_all=matfile(fullfile(Path,'data','single_data',[animal '_all_data.mat']));


            select_id_3{curr_group}=(strcmp([data_all.task_name],'stim_wheel_right_stage1')|...
                strcmp([data_all.task_name],'stim_wheel_right_stage2'))&...
                ~cellfun(@isempty ,data_all.wf_lcr_passive)&...
                ~cellfun(@isempty ,data_all.wf_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_lcr_passive);

            select_id_3{3-curr_group}=(strcmp([data_all.task_name],'stim_wheel_right_stage1_audio_volume')|...
                strcmp([data_all.task_name],'stim_wheel_right_stage2_audio_volume'))&...
                ~cellfun(@isempty ,data_all.wf_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.wf_lcr_passive) &...
                ~cellfun(@isempty ,data_all.face_hml_passive_audio)&...
                ~cellfun(@isempty ,data_all.face_lcr_passive);


            switch curr_mode
                case 1
                    temp_face_passive=data_all.face_lcr_passive;
                    temp_wf_passive_0=data_all.wf_lcr_passive;
                case 2
                    temp_face_passive=data_all.face_hml_passive_audio;
                    temp_wf_passive_0=data_all.wf_hml_passive_audio;
            end

            % pupil

            % temp_idd=find(~cellfun(@isempty,temp_face_passive,'UniformOutput',true));
            % pupil_idx=ismember(1:numel(temp_face_passive),temp_idd(cellfun(@(x)  x.validation.pupil, ...
            %     temp_face_passive(temp_idd),'UniformOutput',true)))';


            temp_behavior=data_all.behavior_task;
            temp_p_val=cellfun(@(x) ...
                arrayfun(@(id)  temp_behavior{id}.rxn_l_p(1)<0.05, find(x),'UniformOutput',true),...
                select_id_3,'UniformOutput',false);



            temp_pupil_size=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx)  feval(@(d) cat(3,d{:}), ...
                cellfun(@(c) nanmean(c,1), xx.pupil_data.diameterZ_filt_sav,'uni',false)),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x ),select_id_3,'UniformOutput',false),'UniformOutput',false);



            pupil_data=cat(4,nan(size(temp_pupil_size{1},1), size(temp_pupil_size{1},2),size(temp_pupil_size{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_pupil_size{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_pupil_size{2}(:,:,:,1:min(5,end)));

            pupil_3stage={temp_pupil_size{1}(:,:,:, temp_p_val{1}==0),...
                temp_pupil_size{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_pupil_size{2}(:,:,:, find(temp_p_val{2}==1,2,'last'))};


            %pupil center

            temp_pupil_center=cellfun(@(aa)  feval(@(a) cat(5,a{:}), cellfun(@(xx)  feval(@(d) cat(4,d{:}), ...
                cellfun(@(c) nanmean(c,1), xx.pupil_data.center_filt_sav,'uni',false)),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x),select_id_3,'UniformOutput',false),'UniformOutput',false);



            pupil_center_data=cat(5,nan(size(temp_pupil_center{1},1), size(temp_pupil_center{1},2),...
                size(temp_pupil_center{1},3),size(temp_pupil_center{1},4),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_pupil_center{1}(:,:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_pupil_center{2}(:,:,:,:,1:min(5,end)));


            pupil_center_3stage={temp_pupil_center{1}(:,:,:,:,temp_p_val{1}==0),...
                temp_pupil_center{1}(:,:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_pupil_center{2}(:,:,:,:,find(temp_p_val{2}==1,2,'last'))};



            % nose
            temp_nose_passive=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx) ...
                permute( nanmean(...
                feval(@(d) cat(5,d{:}), cellfun(@(c) nanmean(c(:,:,3,:),1), xx.face_data.nose_filt_sav,'uni',false)) ,1),[2,5,4,1,3]) ,...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_face_passive(x),select_id_3,'UniformOutput',false),'UniformOutput',false);

            nose_data=cat(4,nan(size(temp_nose_passive{1},1), size(temp_nose_passive{1},2),size(temp_nose_passive{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_nose_passive{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_nose_passive{2}(:,:,:,1:5));

            nose_3stage={temp_nose_passive{1}(:,:,:,temp_p_val{1}==0),...
                temp_nose_passive{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_nose_passive{2}(:,:,:,find(temp_p_val{2}==1,2,'last'))};



            % wf_passive
            temp_wf_passive=cellfun(@(aa)  feval(@(a) cat(4,a{:}), cellfun(@(xx) cat(3,xx.kernels_decoding{:}),...
                aa,'UniformOutput',false)) ,...
                cellfun( @(x) temp_wf_passive_0(x),select_id_3,'UniformOutput',false),'UniformOutput',false);

            wf_passive_3stage={temp_wf_passive{1}(:,:,:, temp_p_val{1}==0),...
                temp_wf_passive{1}(:,:,:, find(temp_p_val{1}==1,2,'last')),...
                temp_wf_passive{2}(:,:,:, find(temp_p_val{2}==1,2,'last'))};


            wf_passive_data=cat(4,nan(size(temp_wf_passive{1},1), size(temp_wf_passive{1},2),size(temp_wf_passive{1},3),...
                5-length(find(temp_p_val{1}==1,5))),...
                temp_wf_passive{1}(:,:,:,find(temp_p_val{1}==1,5,'last')),...
                temp_wf_passive{2}(:,:,:,1:5));




            temp_data_all.nose_passive{curr_animal}=nose_data;
            temp_data_all.pupil_size{curr_animal}=pupil_data;
            temp_data_all.pupil_size_single{curr_animal}=pupil_3stage;
            temp_data_all.pupil_all_trace{curr_animal}=temp_pupil_size;

            temp_data_all.wf_passive{curr_animal}=wf_passive_data;
            temp_data_all.wf_passive_single{curr_animal}=wf_passive_3stage;
            temp_data_all.wf_passive_all_trace{curr_animal}=temp_wf_passive;

            temp_data_all.pupil_center{curr_animal}=pupil_center_data;
            temp_data_all.pupil_center_single{curr_animal}=pupil_center_3stage;
            temp_data_all.pupil_center_all_trace{curr_animal}=temp_pupil_center;

            temp_data_all.nose{curr_animal}=nose_data;
            temp_data_all.nose_single{curr_animal}=nose_3stage;
            temp_data_all.nose_all_trace{curr_animal}=temp_nose_passive;


            clearvars('-except',preload_vars{:});
        end

        all_data.([groups_name{curr_group} '_' modes_name{curr_mode}])=temp_data_all;



    end

end


passive_id={[3 1 2],[2 1 3]};
figure('Position',[50 50 1000 800])
mainfig=tiledlayout( 3,3, ...
    'TileSpacing', 'tight', 'Padding', 'none');

plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 1;  % 明确放在主 layout 的第 1 个 tile

% pupil size  trace
colors={[0 0 1;0 0 0],[1 0 0;0 0 0]}
for curr_group=1:2
    for curr_mode=curr_group

        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_size_single{:}));
        temp_pupil_trace=cellfun(@(x) permute(diff(x,1,2),[2,3,4,1]),pupil_data_single,'UniformOutput',false);
        temp_pupil_trace_2=cellfun(@(x) cat(2,x(:,passive_id{curr_mode}(1),:),nanmean(x(:,passive_id{curr_mode}(2:3),:),2)),...
            temp_pupil_trace,'UniformOutput',false);
        % figure('Position',[50 50 400 150]);
        % tiledlayout(1,2)
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_trace_2{curr_state},3),...
                nanstd(temp_pupil_trace_2{curr_state},0,3)./sqrt(size(temp_pupil_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([-0.01 0.03])

            xline(0)
            axis off
        end


    end
end

hold on
line([-0.2 -0.2],[-0.01 0],'Color','k');
line([-0.2 0],[-0.01 -0.01],'Color','k');
text(-0.32, -0.01, '0.01 arb','Rotation',90)
text(-0.2, -0.015, '0.2 s')

% pupil size across days
plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor"  );
plot_fig.Layout.Tile = 2;  % 明确放在主 layout 的第 1 个 tile

colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2
        % pupil
        pupil_data=cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_size{:});
        temp_avg_pupil_data= cat(3,pupil_data(:,:,passive_id{curr_mode}(1),:),nanmean(pupil_data(:,:,passive_id{curr_mode}(2:3),:),3))
        temp_pupil_AUC_mean=feval(@(a)...
            permute(nanmean(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),1),...
            [4,3,2,1]),diff(temp_avg_pupil_data,1,2));

        temp_pupil_AUC_error=feval(@(a)...
            permute(nanstd(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),0,1)./sqrt(size(pupil_data,1)),...
            [4,3,2,1]),diff(temp_avg_pupil_data,1,2));


        nexttile(plot_fig)
        hold on
        ap.errorfill(1:5,temp_pupil_AUC_mean(1:5,:),temp_pupil_AUC_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_pupil_AUC_mean(6:10,:),temp_pupil_AUC_error(6:10,:),[colors_temp{curr_group};0 0 0])

        xlim([1 10])
        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})
        ylim([-0.1 0.4])
        ylabel('Pupil diameter (arb)')

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end


% wf vs pupil size
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 3;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)
for curr_group=1:2
    for curr_mode=curr_group

        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_all_trace{:}));

        temp_pupil_AUC_single=cellfun(@(x) feval(@(a)...
            permute(trapz(1:sum(face_time>0&face_time<1),a(:,face_time>0&face_time<1,:,:),2),...
            [4,3,1,2]),diff(x,1,2)),pupil_data_single,'UniformOutput',false)



        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}))
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,kernels_period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_pupil_AUC_single,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Pupil diameter')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(-0.1, 0.7, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end

set(gca,'Color','none')
axis square
drawnow


plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 4;  % 明确放在主 layout 的第 1 个 tile

% pupil saccade  trace
passive_id={[3 1 2],[2 1 3]};
colors={[0 0 1;0 0 0;1 0 0],[1 0 0;0 0 0;1 0 0]}
for curr_group=1:2
    for curr_mode=curr_group




        pupil_center_single=feval(@(a) arrayfun(@(id) cat(5,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center_single{:}));

        % temp_pupil_center_trace=cellfun(@(x) permute( vecnorm((x-x(:,1,:,:,:)),2,3),[2,3,4,5,1]),pupil_center_single,'UniformOutput',false);
        temp_pupil_center_trace=cellfun(@(x) permute( vecnorm(diff(x,1,2),2,3),[2,4,5,3,1]),pupil_center_single,'UniformOutput',false);

        temp_pupil_center_trace_2=cellfun(@(x) cat(2,x(:,passive_id{curr_mode}(1),:),nanmean(x(:,passive_id{curr_mode}(2:3),:),2)),...
            temp_pupil_center_trace,'UniformOutput',false);


        % figure('Position',[50 50 400 150]);
        % tiledlayout(1,2)
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_center_trace_2{curr_state},3),...
                nanstd(temp_pupil_center_trace_2{curr_state},0,3)./sqrt(size(temp_pupil_center_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([0 0.03])
            xline(0)
            axis off
        end


    end
end
hold on
line([-0.2 -0.2],[0 0.01],'Color','k');
line([-0.2 0],[0 0],'Color','k');
text(-0.42, 0.00, '0.01 arb','Rotation',90)
text(-0.2, -0.005, '0.2 s')



plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor" );
plot_fig.Layout.Tile = 5;  % 明确放在主 layout 的第 1 个 tile

% pupil saccade
passive_id={[3 1 2],[2 1 3]};
% fig1=figure;
% tiledlayout(fig1,2,2,"TileIndexing","columnmajor")
colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2

        % pupil center
        pupil_center_data=cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center{:});

        temp_pupil_center=permute(cat(4,vecnorm(diff(pupil_center_data(:,:,:,passive_id{curr_mode}(1),:),1,2),2,3),...
            nanmean(vecnorm(diff(pupil_center_data(:,:,:,passive_id{curr_mode}(2:3),:),1,2),2,3),4)), [2,4,5,1,3]);


        temp_pupil_center_mean=permute(nanmean(trapz(temp_pupil_center(face_time>0&face_time<1,:,:,:),1),4),[3,2,1]);
        temp_pupil_center_error=permute(nanstd(trapz(temp_pupil_center(face_time>0&face_time<1,:,:,:),1),0,4)./...
            sqrt(size(temp_pupil_center,4)),[3,2,1]);

        nexttile(plot_fig)
        ap.errorfill(1:5,temp_pupil_center_mean(1:5,:),temp_pupil_center_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_pupil_center_mean(6:10,:),temp_pupil_center_error(6:10,:),[colors_temp{curr_group};0 0 0])
        ylim([0 1])
        xlim([1 10])
        ylabel('Pupil saccade (arb)')

        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end

% wf vs pupil saccade
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 6;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)

for curr_group=1:2
    for curr_mode=curr_group

        pupil_center_single=feval(@(a) arrayfun(@(id) cat(5,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).pupil_center_all_trace{:}));

        temp_pupil_center=cellfun(@(x) permute(vecnorm(diff(x,1,2),2,3),[5,4,2,1,3]),...
            pupil_center_single,'UniformOutput',false);

        temp_pupil_center_mean=cellfun(@(x) trapz(x(:,:,face_time>0&face_time<1),3),temp_pupil_center,'UniformOutput',false);

        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}));
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,kernels_period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_pupil_center_mean,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Pupil saccade')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(0, 1, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end

set(gca,'Color','none')
axis square
drawnow


% nose move  trace

plot_fig=tiledlayout(mainfig,2, 2);
plot_fig.Layout.Tile = 7;  % 明确放在主 layout 的第 1 个 tile
colors={[0 0 1;0 0 0],[1 0 0;0 0 0]}
for curr_group=1:2
    for curr_mode=curr_group
        pupil_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:3,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose_single{:}));
        temp_pupil_trace=cellfun(@(x) permute(vecnorm(diff(x,1,1),2,3),[2,1,4,3]),pupil_data_single,'UniformOutput',false);
        temp_pupil_trace_2=cellfun(@(x) cat(1,x(passive_id{curr_mode}(1),:,:),nanmean(x(passive_id{curr_mode}(2:3),:,:),1)),...
            temp_pupil_trace,'UniformOutput',false);
        for curr_state=1:2
            nexttile(plot_fig)
            hold on
            ap.errorfill(face_time(2:end)',nanmean(temp_pupil_trace_2{curr_state},3)',...
                nanstd(temp_pupil_trace_2{curr_state},0,3)'./sqrt(size(temp_pupil_trace_2{curr_state},3)),colors{curr_group})
            xlim([-0.2 1])
            ylim([-0.5 1])
            xline(0)
            axis off
        end
    end
end

hold on
line([-0.2 -0.2],[-0.5 0],'Color','k');
line([-0.2 0],[-0.5 -0.5],'Color','k');
text(-0.42,  -0.4, '0.5 arb' ,'Rotation',90)
text(-0.1, -0.6, '0.2 s')

% nose move across days
passive_id={[3 1 2],[2 1 3]};

plot_fig=tiledlayout(mainfig,2, 2,"TileIndexing","columnmajor" );
plot_fig.Layout.Tile = 8;  % 明确放在主 layout 的第 1 个 tile

colors_temp={[84 130 53]/255,[112  48 160]/255}
for curr_group=1:2
    for curr_mode=1:2
        % pupil
        nose_data=cat(5, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose{:});
        temp_avg_nose_data= cat(2,vecnorm(diff(nose_data(:,passive_id{curr_mode}(1),:,:,:),1,1),2,3),...
            nanmean(vecnorm(diff(nose_data(:,passive_id{curr_mode}(2:3),:,:,:),1,1),2,3),2));
        temp_nose_AUC_mean= permute(nanmean(trapz(1:sum(face_time>0&face_time<1),temp_avg_nose_data(face_time>0&face_time<1,:,:,:,:),1),5),...
            [4,2,3,1]);

        temp_nose_AUC_error=...
            permute(nanstd(trapz(1:sum(face_time>0&face_time<1),temp_avg_nose_data(face_time>0&face_time<1,:,:,:,:),1),0,5)./sqrt(size(nose_data,5)),...
            [4,2,3,1]);


        nexttile(plot_fig)
        hold on
        ap.errorfill(1:5,temp_nose_AUC_mean(1:5,:),temp_nose_AUC_error(1:5,:),[colors_temp{curr_group};0 0 0])
        ap.errorfill(6:10,temp_nose_AUC_mean(6:10,:),temp_nose_AUC_error(6:10,:),[colors_temp{curr_group};0 0 0])

        xlim([1 10])
        xticks([1 5 6 10])
        xticklabels({'-5','-1','0','4'})
        ylim([0 15])
        ylabel('Nose move (arb)')

        xlabel('Day (s)')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        set(gca,'Color','none')
        drawnow

    end
end

% nose move vs wf
plot_fig=tiledlayout(mainfig,1,1);
plot_fig.Layout.Tile = 9;  % 明确放在主 layout 的第 1 个 tile
nexttile(plot_fig)

for curr_group=1:2
    for curr_mode=curr_group

        nose_move=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).nose_all_trace{:}));

        nose_move_trace=cellfun(@(x) permute(vecnorm(diff(x,1,1),2,3),[4,2,1,3]),...
            nose_move,'UniformOutput',false);

        nose_move_trace_mean=cellfun(@(x) trapz(x(:,:,face_time>0&face_time<1),3),nose_move_trace,'UniformOutput',false);

        wf_data_single=feval(@(a) arrayfun(@(id) cat(4,a{:,id}),1:2,'UniformOutput',false),...
            cat(1, all_data.([groups_name{curr_group} '_' modes_name{curr_mode}]).wf_passive_all_trace{:}));
        all_passive_image_pre=cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),...
            x),wf_data_single,'UniformOutput',false);
        temp_wf_passive_plot_tace= cellfun(@(x) ds.make_each_roi(x, length(t_kernels),roi1),all_passive_image_pre,'UniformOutput',false)
        temp_wf_passive_max=cellfun(@(x) permute(max(x(1,kernels_period,:,:) ,[],2),[4,3,2,1])  ,temp_wf_passive_plot_tace,'UniformOutput',false)

        colors{1}={[0 0 1],[0 0 1]}
        colors{2}={[1 0 0],[1 0 0]}

        hold on
        temp_x=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            nose_move_trace_mean,'UniformOutput',false);
        temp_y=cellfun(@(a)  a(:,passive_id{curr_mode}(1)),...
            temp_wf_passive_max,'UniformOutput',false);

        cellfun(@(x,y,z) plot(x,y,...
            'LineStyle','none','Marker','.','MarkerSize',10,'Color',z),temp_x,temp_y,colors{curr_group},'uni',false)
        % xlim([-0.1 0.6])
        % ylim([0 0.00035])
        ylabel('mPFC \Delta F/F_0')
        xlabel('Nose move')
        % title([groups_name{curr_group} '\_' modes_name{curr_mode}])
        [R(curr_group,1), p(curr_group,1)] = corr(cat(1,temp_x{:}), cat(1,temp_y{:}));

        temp_line = polyfit(cat(1,temp_x{:}), cat(1,temp_y{:}), 1);
        x_fit_task = linspace(0, 20, 2);
        y_fit_task = polyval(temp_line, x_fit_task);
        plot(x_fit_task, y_fit_task, '-', 'LineWidth', 2,'Color',colors{curr_group}{1});



    end
end
set(gca,'Color','none')

axis square
drawnow


exportgraphics(gcf, fullfile(Path,...
    'submission_3_NatureCommunications\revisions\revision_figures\eps\Fig_s14.eps'), ...
    'ContentType','vector');

 clearvars('-except',main_preload_vars{:});
