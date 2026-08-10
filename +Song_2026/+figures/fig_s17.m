%% projection distribution
main_preload_vars = who;

allenAtlasPath = fileparts(which('template_volume_10um.npy'));
saveLocation=fullfile(allenAtlasPath,'temp_connect')
 % saveLocation = 'C:\Users\dsong\Documents\temp_connect'; % where to save the data downloaded from the Allen Connectivity dataset 
% allenAtlasPath =  'C:\Users\dsong\Documents\GitHub\osfstorage-archive'; % download from: https://figshare.com/articles/dataset/Modified_Allen_CCF_2017_for_cortex-lab_allenCCF/25365829 
fileName = ''; % leave empty to recompute each time (e.g. load the Allen raw data and sumnmarize it into one matrix), 
 
all_inputRegions={{'VIS'}, {'AUD'}}
corlors={'B','R'}
temp_img_data=cell(2,1)
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
[projectionMatrix_array, projectionMatrixCoordinates_ARA,Boundarys,brain_boudarys,temp_img_data{curr_fig}]=...
bsv.plotConnectivity(experimentImgs, allenAtlasPath, outputRegions, numberOfSlices, numberOfPixels, plane, regionOnly, smoothing, colorLimits, color);

exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s17_' num2str(curr_fig)  '.eps']), ...
    'ContentType','vector');
end
 



clearvars('-except',main_preload_vars{:});

 

figure('Position',[50 50 1400 600]);
tiledlayout(3,10,'TileSpacing','none')
max_y=1033;
max_x=617;
for curr_i=1:10
    AX=nexttile

    % temp_image=  ds.colormap_overlay(temp_img_data{1}{curr_i}'./max(cat(3,temp_img_data{1}{:}),[],'all'),...
    %     temp_img_data{2}{curr_i}'./max(cat(3,temp_img_data{2}{:}),[],'all'),'B','R');
    % image(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_image)


    imagesc(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_img_data{1}{curr_i}'./max(cat(3,temp_img_data{1}{:}),[],'all'))
    colormap(AX,ap.colormap('WB'))
    clim([0 1])


    hold on
    axis image off
    plot(Boundarys{curr_i}(1,:),Boundarys{curr_i}(2,:),'Color','k','LineWidth',2)
    plot(brain_boudarys{curr_i}(:,2),brain_boudarys{curr_i}(:,1),'Color','k','LineWidth',2)
    % plor
    set(gca, 'YDir', 'reverse');
    tempx=xlim;
    xlim([sum(tempx)/2-max_x/2 sum(tempx)/2+max_x/2]);
    tempy=ylim;
    ylim([sum(tempy)/2-max_y/2 sum(tempy)/2+max_y/2]);

end
% colorbar
for curr_i=1:10
   AX= nexttile
    % 
    % temp_image=  ds.colormap_overlay(temp_img_data{1}{curr_i}'./max(cat(3,temp_img_data{1}{:}),[],'all'),...
    %     temp_img_data{2}{curr_i}'./max(cat(3,temp_img_data{2}{:}),[],'all'),'B','R');
    % image(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_image)


    imagesc(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_img_data{2}{curr_i}'./max(cat(3,temp_img_data{2}{:}),[],'all'))
    colormap(AX,ap.colormap('WR'))
    clim([0 1])


    hold on
    axis image off
    plot(Boundarys{curr_i}(1,:),Boundarys{curr_i}(2,:),'Color','k','LineWidth',2)
    plot(brain_boudarys{curr_i}(:,2),brain_boudarys{curr_i}(:,1),'Color','k','LineWidth',2)
    % plor
    set(gca, 'YDir', 'reverse');
    tempx=xlim;
    xlim([sum(tempx)/2-max_x/2 sum(tempx)/2+max_x/2]);
    tempy=ylim;
    ylim([sum(tempy)/2-max_y/2 sum(tempy)/2+max_y/2]);

end
for curr_i=1:10
    nexttile

    temp_image=  ds.colormap_overlay(temp_img_data{1}{curr_i}'./max(cat(3,temp_img_data{1}{:}),[],'all'),...
        temp_img_data{2}{curr_i}'./max(cat(3,temp_img_data{2}{:}),[],'all'),'B','R');
    image(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_image)


    % imagesc(projectionMatrixCoordinates_ARA{curr_i}{1},projectionMatrixCoordinates_ARA{curr_i}{2},temp_img_data{1}{curr_i}'./max(cat(3,temp_img_data{1}{:}),[],'all'))
    % colormap(ap.colormap('WB'))
    % clim([0 1])


    hold on
    axis image off
    plot(Boundarys{curr_i}(1,:),Boundarys{curr_i}(2,:),'Color','k','LineWidth',2)
    plot(brain_boudarys{curr_i}(:,2),brain_boudarys{curr_i}(:,1),'Color','k','LineWidth',2)
    % plor
    set(gca, 'YDir', 'reverse');
    tempx=xlim;
    xlim([sum(tempx)/2-max_x/2 sum(tempx)/2+max_x/2]);
    tempy=ylim;
    ylim([sum(tempy)/2-max_y/2 sum(tempy)/2+max_y/2]);

end

exportgraphics(gcf, fullfile(Path,['figures\eps\Fig s17'  '.eps']), ...
    'ContentType','vector');