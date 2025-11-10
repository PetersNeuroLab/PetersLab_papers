
%%  Neuro street view example script 

saveLocation = 'C:\Users\dsong\Documents\temp_connect'; % where to save the data downloaded from the Allen Connectivity dataset 
allenAtlasPath =  'C:\Users\dsong\Documents\GitHub\osfstorage-archive'; % download from: https://figshare.com/articles/dataset/Modified_Allen_CCF_2017_for_cortex-lab_allenCCF/25365829 
fileName = ''; % leave empty to recompute each time (e.g. load the Allen raw data and sumnmarize it into one matrix), 
 
all_inputRegions={{'VIS'}, {'AUD'}}
corlors={'B','R'}
A5=cell(2,1)
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

[A1,A2,A3,A4,A5{curr_fig}]=bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color);

% exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s5_' num2str(curr_fig)  '.eps']), ...
%     'ContentType','vector');

end


 % Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';


%%
 Path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Song_2025';

group_color={'B','R'}
img_data=cell(2,1)

 length_maxY=max(cellfun(@(a) max(a)-min(a), cellfun(@(x)  x(:,1) ,A4,'UniformOutput',false),'UniformOutput',true));
 length_maxX=max(cellfun(@(a) max(a)-min(a), cellfun(@(x)  x(:,2) ,A4,'UniformOutput',false),'UniformOutput',true));


for curr_group=1:2
    figure('Position',[50 50 1400 200],'Color','w')
    t=tiledlayout(1,10,'Padding','compact','TileSpacing','loose')
    for curr_slice=1:10

        % ax=figure('Color','w','Position',[50 50 1000 1000]);

        ax1=nexttile(t);      % 激活第一个 tile
        hold on
        imagesc(A2{curr_slice}{1},A2{curr_slice}{2},A5{curr_group}{curr_slice}');

         clim([0 1])
        % clim([0 max(vertcat(A5{1}{curr_slice},A5{2}{curr_slice}),[],'all')])
        % clim([0 max(A5{curr_group}{curr_slice},[],'all')])
        axis equal
        axis square
        axis image off
        colormap(ap.colormap(['W' group_color{curr_group}],[],0.6))
        plot(A3{curr_slice}(1,:), ...
            A3{curr_slice}(2,:), ...
            'Color', [0 0 0], 'LineWidth', 2);

        plot(A4{curr_slice}(:,2), ...
            A4{curr_slice}(:,1), ...
            'Color', [0 0 0], 'LineWidth', 2);

        xlim([ max(A4{curr_slice}(:,2))-length_maxX   max(A4{curr_slice}(:,2))])
        ylim([max(A4{curr_slice}(:,1))-length_maxY    max(A4{curr_slice}(:,1))])
        set(gca, 'YDir', 'reverse');


    end
        cb = colorbar;

    % exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s5_' num2str(curr_group)  '.eps']), ...
    %     'ContentType','vector');

    frame1 = getframe(gcf);
    img_data{curr_group} =im2double( frame1.cdata);
end

result = min(img_data{1}, img_data{2});

figure
imshow(result);

exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s5_merge'  '.eps']), ...
    'ContentType','vector');


img_bar=cell(2,1);
vals = linspace(0, 1, 200)';

for curr_bar=1:2
    figure('Position',[50 50 10 500],'Color','w');

    % 绘制伪 colorbar
    switch curr_bar
        case 1
            imagesc([0 1], [0 1], vals);
        case 2
            imagesc([0 1], [1 0], vals);

    end
    colormap(ap.colormap(['W' group_color{curr_bar}],[],0.6));
    axis off; axis xy;
    temp_frame = getframe(gcf);
    img_bar{curr_bar} =im2double( temp_frame.cdata);

end
result = min(img_bar{1}, img_bar{2});

figure
imshow(result);
exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s5_merge_bar'  '.eps']), ...
    'ContentType','vector');

