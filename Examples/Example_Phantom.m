%% Load example phantom data
% This dataset is obtained from a syringe into which pyruvate was injected.
% To mimic in-vivo effects, the shim was set to the scanner default, i.e.
% good $B_0$ homogeneity was not ensured. 
% 
% Owing to the position of the extension line through which pyruvate was
% injected, it is possible to observe pyruvate to the right of the syringe
% body on some (but not all) frames owing to the shifting FOV caused by a
% residual ~50 Hz off-resonance shift. 
%
% The below .mat file contains:
%
%   imagesIn -- a 1x1 cell, containing 4D hyperpolarized {1} pyruvate images 
%   in the order (xyzt). 
%   
%   carbonParams -- some acquisition parameters relating to the EPI data 
%   
%   H1Full -- a high-ish resolution anatomical image 
%   H1Params -- some acquisition parameters for the anatomical image 
%   H1SameRes -- anatomical image at the same resolution as the HP images 
%   
 


% Phantom dataset: 
load ./data/PhantomData_1.mat
peSignConvention=1; 
%                    Which direction the first shot in the image traversed 
%                    k-space; as the scan was started prior to the
%                    pyruvate injection (and proceded with alternating ±
%                    PE gradient signs) it is not well defined when the
%                    first bolus arrives -- this parameter encodes what the
%                    first frame's polarity was. 
%

% In vivo datasets: 
% load ./data/ExampleData_1.mat
% peSignConvention=1; 
%
% Additional datasets: 
% load ./data/ExampleData_2.mat
% peSignConvention=0; 
%
% load ./data/ExampleData_3.mat
% peSignConvention=0;



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
viewimage(cat_images_5d(permute(imagesIn{1},[1 2 3 5 4])),'title','Original Image', 'zoom', 2)
viewimage(cat_images_5d(permute(out{1},[1 2 3 5 4])),'title','Corrected Image', 'zoom', 2)

%% Create an intensity mask based on the anatomical 1H images
% Rather than manually segmenting the syringe, do it based on thresholding

mask=zeros(size(H1SameRes)); 
protonThresh=0.5;
for ix=1:size(mask,3) 
    mask(:,:,ix)=im2bw(uint8(abs(H1SameRes(:,:,ix))),protonThresh); 
end

viewimage(cat_images(mask)); 
%% Compare with the HP images by masking based on their intensity
% This is to provide a quantitative evaluation of the method. 

HPGamma=1;%Optionally do a gamma transform to ameliorate the surface coil 
% profile and the flow artefact visible while the injection proceeds in the
% middle of the image. 
HPThresh=0.5; %Threshold for HP images

HPMaskOriginal=(im2bw(imadjust(uint16(abs(cat_images(sum(imagesIn{1},4)))),[],[],HPGamma),HPThresh)); 
HPMaskCor=(im2bw(imadjust(uint16(abs(cat_images(sum(out{1},4)))),[],[],HPGamma),HPThresh)); 
viewimage(cat_images(HPMaskOriginal))
viewimage(cat_images(HPMaskCor))
%% Compute pyruvate-to-proton Jaccard index
%Jaccard index is defined as the intersection / union of sets
disp('Pyruvate: Corrected and original Jaccard indicies compared to proton'); 
Jaccard_Index_Original = sum(HPMaskOriginal(:) & mask(:)) / sum(HPMaskOriginal(:) | mask(:))
Jaccard_Index_Corrected = sum(HPMaskCor(:) &  mask(:)) / sum(HPMaskCor(:) | mask(:))

%% Do the same thing, comparing even & odd echoes to themselves (temporally summed over time): 

HPMaskOddOriginal=(im2bw(imadjust(uint16(abs(cat_images(sum(imagesIn{1}(:,:,:,1:2:end),4)))),[],[],HPGamma),HPThresh)); 
HPMaskEvenOriginal=(im2bw(imadjust(uint16(abs(cat_images(sum(imagesIn{1}(:,:,:,2:2:end),4)))),[],[],HPGamma),HPThresh)); 
disp('Original Jaccard self-index:')
Jaccard_Index_Self_Original = sum(HPMaskOddOriginal(:) & HPMaskEvenOriginal(:)) / sum( HPMaskOddOriginal(:) | HPMaskEvenOriginal(:)) 


HPMaskOddCorrected=(im2bw(imadjust(uint16(abs(cat_images(sum(out{1}(:,:,:,1:2:end),4)))),[],[],HPGamma),HPThresh)); 
HPMaskEvenCorrected=(im2bw(imadjust(uint16(abs(cat_images(sum(out{1}(:,:,:,2:2:end),4)))),[],[],HPGamma),HPThresh)); 
disp('Corrected Jaccard self-index:'); 
Jaccard_Index_Self_Corrected = sum(HPMaskEvenCorrected(:) & HPMaskOddCorrected(:) ) / sum(HPMaskEvenCorrected(:) | HPMaskOddCorrected(:) )