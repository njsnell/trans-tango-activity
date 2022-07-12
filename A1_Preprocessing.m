% Script 1: Preprocessing
% 
% Requires NoRMCorre and CaImAn matlab packages to run


%% Convert multiple tif stacks into mat file

% Adjust these parameters:
folder_to_analyze = '';         % folder with all imaging .tif files
n_stimuli = 3;
n_trials_per_stimulus = 3;
n_frames_per_trial = 480;   
n_planes = 12;
ly = 256;                       % length in y dimension (# pixels)
lx = 256;                       % length in x dimension (# pixels)
downsample_factor = 2;
use_nonrigid = true;
% Don't change anything below!!!

% compute more parameters
lx = lx/downsample_factor; ly = ly/downsample_factor;   % downsample
nTrials = n_trials_per_stimulus * n_stimuli;
nFrames = n_frames_per_trial / n_planes;
n_frames_total_per_stimulus = n_frames_per_trial*n_trials_per_stimulus;
Frames_per_Plane = n_frames_per_trial * n_trials_per_stimulus * n_stimuli / n_planes;

% get list of tif files
addpath(folder_to_analyze);
folderFilesAll = dir(folder_to_analyze);
firstfile = 3;
if strcmp(folderFilesAll(firstfile).name,'.DS_Store')
    firstfile = 4;
end
folderFilesAll(firstfile).name
folderFiles = folderFilesAll(firstfile:firstfile+nTrials-1);

for i = 1:n_trials_per_stimulus 
    for j = 1:n_stimuli 
        fileIndex = (j-1)*n_trials_per_stimulus + i;
        stackFileName = folderFiles(fileIndex).name;
        path_names{i,j} = [folder_to_analyze stackFileName];
    end
end

for j = 1:n_trials_per_stimulus 
    for i = 1:n_frames_per_trial
        for k = 1:n_stimuli 
            curframe = imread(path_names{j,k},i);
            dsframe = curframe(1:downsample_factor:end,1:downsample_factor:end);
            raw_data(:,:,i,j,k) = dsframe;
        end
    end
end

% Reshape data to have dimensions (x,y,all frames concatenated)
% Add more stimuli as necessary...
stim_1 = reshape(raw_data(:,:,:,:,1), lx, ly, n_frames_total_per_stimulus, []);
stim_2 = reshape(raw_data(:,:,:,:,2), lx, ly, n_frames_total_per_stimulus, []);
stim_3 = reshape(raw_data(:,:,:,:,3), lx, ly, n_frames_total_per_stimulus, []);
all_stim = ...
    cat(3, stim_1, stim_2, stim_3); 

% Reshape data to have dimensions: (x,y,frames per plane,planes)
planes = zeros(ly,lx,Frames_per_Plane,n_planes);
for i = 1:n_planes
    planes(:,:,:,i) = ...
        cat(3, all_stim(:,:,i:n_planes:end));
end

mocoPlanesAll = NaN(n_planes,nTrials,nFrames,lx,ly);    % all motion-corrected planes
for p = 1:n_planes
    rawPlane = squeeze(planes(:,:,:,p));
    fprintf(['Formatting plane ' num2str(p) '...\n']);
    mocoPlaneStruct = Motion_Correction(rawPlane,use_nonrigid);
    mocoPlane = mocoPlaneStruct.M;
    mocoPlane = shiftdim(mocoPlane,2);
    for t = 1:nTrials
        frameRange = (t-1)*nFrames + [1:nFrames];
        mocoPlanesAll(p,t,:,:,:) = mocoPlane(frameRange,:,:);
    end
end

%% Save
save([folder_to_analyze 'mocoPlanesAll.mat'],'-v7.3','mocoPlanesAll');
close all;

%% Preprocess fly data
stimOn_frameList = [];      % Array of indices of frame of stimulus onset for each trial
stimOff_frameList = [];     % Array of indices of frame of stimulus offset for each trial    
OnOff = 0;                          % 1 if analyzing ON and OFF; 0 otherwise
ly = 128;
lx = 128;
n_planes = 12;
nTrials = 9; 
nFrames_Output = 20;                % use 20 for On only, 25 for On/Off

if OnOff
    data = AlignToStimulus_OnOff(mocoPlanesAll,stimOn_frameList,stimOff_frameList);
else
    data = AlignToStimulus(mocoPlanesAll,stimOn_frameList);
end

%% Extract ROIs using CaImAn's CNMF algorithm
% Adapted from CaImAn package

Y = squeeze(data);
if ndims(Y) == 4
    [d1,d2,d3,T] = size(Y);                            % dimensions of dataset
else
    [d1,d2,T] = size(Y);
    d3 = 1;
end
d = d1*d2*d3;                                          % total number of pixels
tau = [4 4 2];      % [4,4,2];
k = 200;            % 200;
merge = 0.9;        % 0.9;
p = 0;                                            % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'d3',d3,...                  % dimensions of datasets
    'search_method','dilate',...                 % search locations when updating spatial components
    'maxIter',15,...                             % number of NMF iterations during initialization
    'deconv_method','constrained_foopsi',...     % activity deconvolution method
    'temporal_iter',2,...                        % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge,...                            % merging threshold
    'rem_prct',10,...
    'gSig',tau,'nb',1 ...
    );

[P,Y] = preprocess_data(Y,p);
Cn = correlation_image_3D(Y); % for large datasets change with reshape(P.sn,d1,d2,d3), %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)

[Ain,Cin,bin,fin,center] = initialize_components(Y,k,tau,options,P);  % initialize
ff = find(sum(Ain)<1e-3*mean(sum(Ain)));   % remove very small components
Ain(:,ff) = [];
Cin(ff,:) = [];
center(ff,:) = [];

Yr = reshape(Y,d,T);

[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

P.p = 0;
[C,f,P,S,~] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

[A,C,~,~,P,~] = merge_components(Yr,A,b,C,f,P,S,options);

[C_df,~] = extract_DF_F(Yr,A,C,P,options);

SC = A;     % spatial components
TC = C_df;  % temporal components
