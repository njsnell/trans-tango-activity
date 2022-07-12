% Script for automatically analyzing ROI data

% You need to have the following data structures loaded:
% SC_all = 1xNFlies cell array, where each entry is one fly's spatial components (SC) array (obtained from CaImAn) 
% TC_all = 1xNFlies cell array, where each entry is one fly's temporal components (TC) array (obtained from CaImAn) 

%% Input general parameters here
onOff = 0;
offDelays = [];     % Leave as empty array - If looking at both ON and OFF, this value will be computed in next section
nFlies = 8;
nStimTypes = 4;
sigma = 2;
d1 = 128; d2 = 128; d3 = 12;
color1 = [1 0 0];   % in RGB values
color2 = [0 1 0];
color3 = [0 0 1];
colors = [color1; color2; color3];
frames_per_trial = 20;       % number of frames analyzed when converted to ROIs - should be same as size(TC,2)

%% If looking at ON and OFF - Input index of all ON and OFF frames
% To calculate offDelays - number of frames between onset to offset
all_on_frames = ... % [Fly 1 on frames; Fly 2 on frames; etc.]
    [];
all_off_frames = ... % [Fly 1 off frames; Fly 2 off frames; etc.]
    [];
all_off_delays = all_off_frames-all_on_frames;

%% Input ROI tuning parameters here
if onOff
    % Looking at ON and OFF responses
    [~,respType] = MakeCType(2*nStimTypes); % Don't adjust
    stim1_ON = respType{1}; stim1_OFF = respType{2}; % Include an ON and OFF for each stimulus type
    stim2_ON = respType{3}; stim2_OFF = respType{4};
    stim3_ON = respType{5}; stim3_OFF = respType{6};
else
    % Looking at ON responses only
    [~,respType] = MakeCType(nStimTypes); % Don't adjust 
    stim1 = respType{1}; stim2 = respType{2}; stim3 = respType{3}; % Include one for each stimulus type
end

% % Specify ROI tunings here. Use setdiff and intersect functions.
tuning1 = setdiff([stim1 stim2],[stim3 stim4]); % stim 1 or stim 2, and NOT stim 3 or stim 4
tuning2 = setdiff([stim3 stim4],[stim1 stim2]); % stim 3 or stim 4, and NOT stim 1 or stim 2
tuning3 = intersect([stim1 stim2],[stim3 stim4]); % stim 1 or stim 2, AND stim 3 or stim 4

% Put all tunings into a cell array.
tuning = {tuning1,tuning2,tuning3};

%% Count ROIs
clCounts = zeros(nFlies,length(tuning));
for fly = 1:nFlies
    TC = TC_all{fly};
    if onOff
        offDelays = all_off_delays(fly,:);
    end
    [traces,~] = NormalizeTraces(TC,nStimTypes,frames_per_trial);
    c_out = SeparateTracesBySelectivityGroups(traces,sigma,nStimTypes,frames_per_trial,offDelays);
    for c = 1:length(tuning)
        idx_c = ismember(c_out,tuning{c});
        clCounts(fly,c) = sum(idx_c);
    end
end
figure; b = bar(squeeze(clCounts(:,:,1)),'stacked'); 
b(1).FaceColor = color1; b(2).FaceColor = color2; b(3).FaceColor = color3;
yticks(20:20:140);
xlabel('Fly #'); ylabel('# ROIs per tuning');
f = gcf; a = gca;
a.FontSize = 20; a.Box = 'off'; f.Color = 'w';

%% Compute stats: intersection over union
close all;
if onOff
    matrixSize = 2*nStimTypes;
    r_perm = [1:2:nStimTypes*2, 2:2:nStimTypes*2];
else
    matrixSize = nStimTypes;  
    r_perm = 1:nStimTypes;
end
all_clIoU = zeros(nFlies,matrixSize,matrixSize);
mean_clIoU = zeros(matrixSize,matrixSize);
for fly = 1:nFlies
    clIoU = zeros(matrixSize,matrixSize);
    TC = TC_all{fly};
    if onOff
        offDelays = all_off_delays(fly,:);
    end
    [traces,~] = NormalizeTraces(TC,nStimTypes,frames_per_trial);
    c_out = SeparateTracesBySelectivityGroups(traces,sigma,nStimTypes,frames_per_trial,offDelays);
    for ix = 1:matrixSize
        for jx = 1:matrixSize
            i = r_perm(ix); j = r_perm(jx);
            tuning_i = intersect(respType{i},respType{j});
            tuning_u = union(respType{i},respType{j});
            c_i = sum(ismember(c_out,tuning_i),'all');
            c_u = sum(ismember(c_out,tuning_u),'all');
            clIoU(ix,jx) = c_i / c_u;
        end
    end
    all_clIoU(fly,:,:) = clIoU;
end
mean_clIoU(:,:) = mean(all_clIoU,1,'omitnan');
fig1 = figure; 
imagesc(mean_clIoU); colormap jet; caxis([0 1]); colorbar;
xticks([]); yticks([]);
fig1.Color = 'w';
fig1.Children(1).Ticks = [0 0.5 1];
fig1.Children(1).FontSize = 22;

%% Make ROI maps for subsections of Z
close all;
nSect = 2;
nClusters = numel(tuning);
allZMaps = NaN(nFlies,nClusters,d1,d2,3,nSect);
for fly = 1:nFlies
    TC = TC_all{fly}; SC = SC_all{fly};
    if onOff
        offDelays = all_off_delays(fly,:);
    end
    [traces,removed] = NormalizeTraces(TC,nStimTypes,frames_per_trial);
    c_out = SeparateTracesBySelectivityGroups(traces,sigma,nStimTypes,frames_per_trial,offDelays);
    SC(:,removed) = [];
    typeClusters = zeros(size(c_out,1),1);
    for ct = 1:length(tuning)
        typeClusters(ismember(c_out,tuning{ct})) = ct;
    end
    for cl = 1:nClusters
        [map1,map2] = AutoselectROIs(SC,typeClusters,cl,nSect,colors);
        allZMaps(fly,cl,:,:,:,1) = map1; 
        allZMaps(fly,cl,:,:,:,2) = map2; 
    end
end

%% Graph ROI maps
save_folder = []; % Leave empty if you don't want to save maps

close all;
for sec = 1:2
    if sec == 1
        type = 'Ant';
    else
        type = 'Post';
    end
    range = (sec-1)*6 + (1:6);
    mapsToGraph = squeeze(allZMaps(:,:,:,:,:,sec)); 
    alpha = 0.9;
    for fly = 1:size(mapsToGraph,1)
        fig2 = figure; fig2.Position = [400 250 400 400];
        curmap = squeeze(mapsToGraph(fly,1,:,:,:));
        imagesc(curmap);
        hold on;
        for cl = 2:nClusters
            curmap = squeeze(mapsToGraph(fly,cl,:,:,:));
            alphamap = alpha .* ceil(max(curmap,[],3));
            imagesc(curmap,'AlphaData',alphamap);        
        end 
        xticks([]); yticks([]);
        if ~isempty(save_folder)
            filenameFig = [save_folder 'ROI Map ' type ' Fly ' num2str(fly)]; 
            savefig(fig2,[filenameFig '.fig']); 
        end
    end
end

%% Manually choose individual ROIs for analysis
% Instructions:
% Select one (and only one) ROI from the map using the crosshairs.
% When selected, the ROI's border will turn white, and calcium trace will
% display to the right.
% Click 'done' in the second window when ROI is selected.
% Then click anywhere in the ROI map to proceed to the next fly.

% Enter parameters:
ROI_type = 1;       % Provide number of ROI type, if looking at multiple ROI types
section = 'a';      % 'a' or 'p' for anterior or posterior sections 

% Don't change below:
close all;
for fly = 1:nFlies
    fprintf(['Analyzing Fly ' num2str(fly) '.\n']);
    TC = TC_all{fly}; SC = SC_all{fly};
    if onOff
        offDelays = all_off_delays(fly,:);
    end
    [traces,removed] = NormalizeTraces(TC,nStimTypes,frames_per_trial);
    c_out = SeparateTracesBySelectivityGroups(traces,sigma,nStimTypes,frames_per_trial,offDelays);
    SC(:,removed) = [];
    typeClusters = zeros(size(c_out,1),1);
    for ct = 1:length(tuning)
        typeClusters(ismember(c_out,tuning{ct})) = ct;
    end
    [out_map,out_trace] = SelectROIs(SC,traces,typeClusters,section,size(tuning,2));
    if ~isempty(out_trace)
        allFlyMaps(fly,ROI_type,:,:,:) = out_map;
        allTraces(fly,ROI_type,:) = out_trace;
    end
end

%% Make mean trace plots
nROItypes = size(allTraces,2);
colors = hsv(nROItypes);
for i = 1:nROItypes
    all_meanT = zeros(nFlies,frames_per_trial*nStimTypes);
    for fly = 1:nFlies
        T = squeeze(allTraces(fly,i,:));
        meanT = ComputeMeanTrace(T,frames_per_trial);
        all_meanT(fly,:) = meanT;
    end
    
    f = figure; f.Position = [1100 345 860 360];
    d_m = mean(all_meanT,1,'omitnan');
    d_e = 1.96*std(all_meanT,0,1,'omitnan')/sqrt(size(all_meanT,1));
    alpha = 0.3; x = 1:size(all_meanT,2); x_vec = [x fliplr(x)];
    hold on;
    patch = fill(x_vec,[(d_m+d_e) fliplr(d_m-d_e)], colors(i,:));
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', alpha);
    plot(x,d_m,'color',colors(i,:),'LineWidth',2);
    a = gca;
    a.YAxisLocation = 'right'; a.XColor = 'w'; a.YColor = 'k';
    xticks([]); yticks([1 4.5]); yticklabels({0 350});
    a.YLim(1) = 0.6; a.YLim(2) = 4.5;
    a.FontSize = 28;
    ylabel(['% \DeltaF/F']);
    set(f,'color','w');   
end

%% Graph/save maps of manually chosen ROIs
colors_manualROIs = []; % Input RGB values of colors to plot manually chosen ROIs in
save_folder = []; % Leave empty if you don't want to save maps

close all;
whiteMap = NaN(nFlies,128,128,3);
for fly = 1:nFlies
    TC = TC_all{fly}; SC = SC_all{fly};
    if onOff
        offDelays = all_off_delays(fly,:);
    end
    [traces,removed] = NormalizeTraces(TC,nStimTypes,frames_per_trial);
    c_out = SeparateTracesBySelectivityGroups(traces,sigma,nStimTypes,frames_per_trial,offDelays);
    SC(:,removed) = [];
    typeClusters = zeros(size(c_out,1),1);
    for ct = 1:length(tuning)
        typeClusters(ismember(c_out,tuning{ct})) = 1;
    end
    map = AutoselectROIs(SC,typeClusters,1,1,[1 1 1]);
    whiteMap(fly,:,:,:) = map;
end
alpha = 1;
for fly = 1:size(allFlyMaps,1)
    % Recolor allFlyMaps
    newmap = zeros(128*128,size(allFlyMaps,2));
    for m = 1:size(allFlyMaps,2)
        map_in = squeeze(allFlyMaps(fly,m,:,:,:));
        map_w = max(map_in,[],3);
        newmap(:,m) = reshape(map_w,128*128,1);
    end
    newmap_color = newmap * colors_manualROIs;
    curmap = reshape(newmap_color,128,128,3);
    % Plot on white BG
    f = figure; f.Position = [100+20*fly 300 400 400];
    bgmap = squeeze(whiteMap(fly,:,:,:));
    image(bgmap);
    hold on;
    alphamap = alpha .* ceil(max(curmap,[],3));
    image(curmap,'AlphaData',alphamap);     
    drawnow;
    xticks([]); yticks([]);
    if ~isempty(save_folder)
        filenameFig = [save_folder 'Manual ROI Examples Fly ' num2str(fly)]; 
        savefig(f,[filenameFig '.fig']); 
    end
end