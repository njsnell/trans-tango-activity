
%% Setup parameters

folder_to_analyze = [];
stimOn_frameList = [];
stimConditions = {'A','B','C','D'};
nTrialsPerStim = 3;  
nPlanes = 12;
xRes = 128;
yRes = 128;

%% Generate heatmaps
% Probably don't change anything below this
load([folder_to_analyze 'mocoPlanesAll.mat']);

imageData = full(mocoPlanesAll);
nStimConditions = length(stimConditions); 

preStimBuffer = 1;
preStimRange = 10; 
postStimBuffer = 1;
postStimRange = 3;
nTrials = nTrialsPerStim * nStimConditions;
preFirstFrames = stimOn_frameList - (preStimBuffer + preStimRange - 1);
preLastFrames = stimOn_frameList - preStimBuffer;
postFirstFrames = stimOn_frameList + postStimBuffer;
postLastFrames = stimOn_frameList + (postStimBuffer + postStimRange - 1);
baselineIntervals = [preFirstFrames' preLastFrames'];
analysisIntervals = [postFirstFrames' postLastFrames'];

bgsub = prctile(imageData,70,'all');    
allBorderPixels = zeros(nPlanes,yRes,xRes);

for p = 1:nPlanes   % compute border pixels
    borderPixels = zeros(yRes,xRes);
    for t = 1:nTrials
        maxOverFrames = squeeze(max(imageData(p,t,:,:,:),[],3));
        borderPixels = borderPixels | (maxOverFrames == 0);
    end
    allBorderPixels(p,:,:) = borderPixels;
end

for t = 1:nTrials
    baselineInterval = baselineIntervals(t,:);
    baselineRange = baselineInterval(1):baselineInterval(2);
    analysisInterval = analysisIntervals(t,:);
    analysisRange = analysisInterval(1):analysisInterval(2);
    for p = 1:nPlanes
        curData = squeeze(imageData(p,t,:,:,:)) - bgsub;     % background subtraction
        curBorderPixels = squeeze(allBorderPixels(p,:,:));      % remove border pixels
        for f = [baselineRange analysisRange]
            temp = squeeze(curData(f,:,:));
            temp(curBorderPixels==1) = 0;
            curData(f,:,:) = temp;
        end

        baselineFrames = curData(baselineRange,:,:); 
        F0 = squeeze(mean(baselineFrames,1));
        F0(F0<0) = 0;
        std0 = squeeze(std(baselineFrames,0,1));
        F = curData(analysisRange,:,:);

        alldF = zeros(postStimRange,yRes,xRes);
        for f = 1:postStimRange
            curF = squeeze(F(f,:,:));
            curdF = 100*(curF - F0) ./ F0;    % dF/F0
            curdF(isinf(curdF)) = 0; curdF(isnan(curdF)) = 0;
            curdF = medfilt2(curdF);        % apply median filter - removes salt and pepper noise
            alldF(f,:,:) = curdF;
        end  
        mean_dF(:,:) = mean(alldF,1);       % median across time
        mean_dF = imgaussfilt(mean_dF,1);
        all_dF_F0(p,t,:,:) = mean_dF;
    end

    curTrialFrames = squeeze(all_dF_F0(:,t,:,:));
    fullHM = squeeze(max(curTrialFrames,[],1));
    trial = mod(t-1,3) + 1;
    stim = stimConditions{ceil(t/3)};
    filename = [stim '-' int2str(trial)];
    savefolder = [folder_to_analyze 'Heatmaps/'];
    if ~exist(savefolder,'dir'); mkdir(savefolder); end
    savefolderImg = [folder_to_analyze 'Images/'];
    if ~exist(savefolderImg,'dir'); mkdir(savefolderImg); end

    imgdata = uint32(fullHM);
    tif_output = Tiff([savefolder filename '.tif'],'w');
    tagstruct.ImageLength = size(imgdata,1);
    tagstruct.ImageWidth = size(imgdata,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    setTag(tif_output,tagstruct);
    write(tif_output,imgdata);

end
