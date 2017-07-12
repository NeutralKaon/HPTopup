%% Load example data
% Each EPI dataset is obtained from a fed healthy rat.
% The below .mat file contains:
%
%   imagesIn -- a 1x3 cell, containing 4D hyperpolarized {1} pyruvate, {2}
%   bicarbonate and {3} lactate images in the order (xyzt). The number of
%   frames provided varies between the datasets, but the matrix size is
%   constant.
%
%   carbonParams -- some acquisition parameters relating to the EPI data
%
%   H1Full -- a high-ish resolution anatomical image
%   H1Params -- some acquisition parameters for the anatomical image
%   H1SameRes -- anatomical image at the same resolution as the HP images
%

% In vivo datasets:
% 1: 
%load ./data/ExampleData_1.mat
%peSignConvention=1;
%HPThresh=[0.25 0.7 0.8];

%2:
%load ./data/ExampleData_2.mat
%peSignConvention=0;
%HPThresh=[0.2 0.7 0.6]; 

%3: 
load ./data/ExampleData_3.mat
peSignConvention=0;
HPThresh=[0.2 0.75 0.6]; %Threshold for HP images

% NB: peSignConvention = ±1 denotes which direction the first shot in the image traversed
%                    k-space; as the scan was started prior to the
%                    pyruvate injection (and proceded with alternating ±
%                    PE gradient signs) it is not well defined when the
%                    first bolus arrives -- this parameter encodes what the
%                    first frame's polarity was.
%
%   HPThresh -- parameter used for estimating the method's performance by
%   thresholding the hyperpolarized images and comparing similarity between
%   odd and even echoes (which, ignoring biology, should be the same)
%
%   HPGamma -- owing to the surface coil profile, this process is aided by
%   the introduction of a gamma transformation prior to thresholding.


%% Create a temporary folder for storing the results in
dn = tempname;
mkdir(dn);
%% Run method
out = RunHPTopup(imagesIn, peSignConvention, dn, carbonParams);
%% Compare images
% See https://github.com/NeutralKaon/VarianTools/tree/master/Misc/ for
% viewimage.m and https://github.com/NeutralKaon/VarianTools/tree/master/Recon/utils
% for cat_images [_5d] .m

%'1' is pyruvate
%'2' is bicarbonate
%'3' is lactate
viewimage(cat_images_5d(permute(imagesIn{2},[1 2 3 5 4])),'title','Original Image', 'zoom', 2)
viewimage(cat_images_5d(permute(out{2},[1 2 3 5 4])),'title','Corrected Image', 'zoom', 2)

%% Simply compare even & odd echoes to themselves (temporally summed over time):
metNames={'Pyruvate','Bicarbonate','Lactate'};
close all; 

for ix = 1:3
    
    HPGamma=1;
    % Optionally do a gamma transform to ameliorate the surface coil
    % profile and the flow artefact visible while the injection proceeds in the
    % middle of the image.
    
    
    fprintf(1,'For metabolite %s:\n', metNames{ix});
    tmp=abs(cat_images(sum(imagesIn{ix}(:,:,:,1:2:end),4))); 
    HPMaskOddOriginal=im2bw(imadjust(uint16(65535*tmp./max(tmp(:))),[],[],HPGamma),HPThresh(ix));
    
    
    tmp=abs(cat_images(sum(imagesIn{ix}(:,:,:,2:2:end),4))); % This is to get around odd im2uint16 behaviour 
    HPMaskEvenOriginal=im2bw(imadjust(uint16(65535*tmp./max(tmp(:))),[],[],HPGamma),HPThresh(ix));
    disp('Original Jaccard self-index:')
    Jaccard_Index_Self_Original = sum(HPMaskOddOriginal(:) & HPMaskEvenOriginal(:)) / sum( HPMaskOddOriginal(:) | HPMaskEvenOriginal(:))
    % Display union of both masks
    viewimage(HPMaskEvenOriginal|HPMaskOddOriginal)
    
    % Do the same thing but for topup
    tmp=abs(cat_images(sum(out{ix}(:,:,:,1:2:end),4))); 
    HPMaskOddCorrected=im2bw(imadjust(uint16(65535*tmp./max(tmp(:))),[],[],HPGamma),HPThresh(ix));
    
    tmp=abs(cat_images(sum(out{ix}(:,:,:,2:2:end),4))); 
    HPMaskEvenCorrected=im2bw(imadjust(uint16(65535*tmp./max(tmp(:))),[],[],HPGamma),HPThresh(ix));
    disp('Corrected Jaccard self-index:');
    Jaccard_Index_Self_Corrected = sum(HPMaskEvenCorrected(:) & HPMaskOddCorrected(:) ) / sum(HPMaskEvenCorrected(:) | HPMaskOddCorrected(:) )
    viewimage(HPMaskEvenCorrected|HPMaskOddCorrected) 
end