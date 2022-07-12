function Mstruct = Motion_Correction(rawPlane,useNonrigid)

Y = single(rawPlane);                                                               % convert to single precision 
T = size(Y,ndims(Y));

% set parameters 
options1 = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'us_fac',50,'init_batch',100,'max_shift',8,'boundary','zero');
options2 = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'us_fac',50,'init_batch',100,'max_shift',64,'boundary','zero');
options3 = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'us_fac',50,'init_batch',100,'max_shift',8,'boundary','zero',...
    'grid_size',[32,32],'mot_uf',4,'max_dev',4);
options4 = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'us_fac',50,'init_batch',100,'max_shift',64,'boundary','zero',...
    'grid_size',[32,32],'mot_uf',4,'max_dev',4);

% perform motion correction
tic; [M1,shifts1,template1,options1] = normcorre_batch(Y,options1); toc;        % Rigid, small shift
tic; [M2,shifts2,template2,options2] = normcorre_batch(Y,options2); toc         % Rigid, big shift
if useNonrigid
    tic; [M3,shifts3,template3,options3] = normcorre_batch(Y,options3); toc     % Nonrigid, small shift
    try
        tic; [M4,shifts4,template4,options4] = normcorre_batch(Y,options4); toc % Nonrigid, big shift
    catch
        M4 = M3;
    end
else
    M3 = M1;
    M4 = M2;
end

[cM1,mM1,vM1] = motion_metrics(M1,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
[cM3,mM3,vM3] = motion_metrics(M3,10);
[cM4,mM4,vM4] = motion_metrics(M4,10);

matrixCells = {M1 M2 M3 M4};
meanArray = [mean(cM1),mean(cM2),mean(cM3),mean(cM4)];
[~,maxI] = max(meanArray);  % Choose motion correction version (rigid/nonrigid, small/big shift) with best average correlation coefficient
M = matrixCells{maxI};

Mstruct.M = M;
