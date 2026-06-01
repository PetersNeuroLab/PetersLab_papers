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
