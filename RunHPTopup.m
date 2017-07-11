function out = RunHPTopup(imagesIn, peSignConvention, dn, params, varargin)
% Calculate and apply distortion correction for 3D hyperpolarized EPI data
% acquired with an alternating phase encoding gradient.
%
% Inputs:
%
% -- imagesIn should be a 1 x nmetabolite-long cell array containing 4D (xyzt)
% images acquired with alternating phase encoding blips. It is assumed that
% each metabolite image has the same matrix size (and an error will result
% if this isn't the case). Data can be complex.
%
% -- if peSignConvention == 0, it is assumed that the PE phase encoding sign
% is [1 -1] alternating between frames; if it is 1, the convention is
% assumed to be reversed (i.e. [-1 1]).
%
% -- dn is the directory name of a working directory in which to store
% distorted/original images, estimated fields and transformation matricies
% (all in nifti/analyse format)
%
% -- params is a struct containing acquisition parameters. These must
% include, at a minimum:
%       params.lro -- readout axis length (in cm)
%       params.lpe -- in-plane phase encoding axis (cm)
%       params.lpe2 -- through plane phase encoding axis (cm)
%       params.etl -- echo train length of the EPI acquisition
%       params.esp -- echo spacing of the EPI acquisition (in s)
%
% -- Optionally, an 'options' string can be provided (as varargin) and is
% passed through as a command line argument to the FSL tool topup. 
%
% FSL is a prerequisite
%
% (c) Jack J. Miller, University of Oxford, 2015-17
% Released under the GPL


useWeights=1;
%% Input checking
setenv('FSLDIR','/usr/local/fsl');


if ~iscell(imagesIn)
    error('Input is not a cell')
end
if ~ismatrix(imagesIn)
    error('Check dimensions of input')
end

sizes=cellfun(@size,imagesIn,'UniformOutput',false);

if length(imagesIn) > 1
    if ~isequal(sizes{:})
        error('Every metabolite image stack must have the same dimensions -- add zeros if required');
    end
end


%% Create equal SNR image
%Assume error constant across metabolites; renorm'd sum:

zfsTest=zeros([sizes{1}(1:3) 2]);
snrweights = zeros(sizes{1}(3),sizes{1}(4),length(imagesIn)); %Image weights


nMets=length(imagesIn);%Number of metabolites

for mx = 1:nMets %Loop over metabolites
    
    
    wt = ones(sizes{1}(1:4));
    if useWeights
        for ix = 1:sizes{1}(4)
            for jx = 1:sizes{1}(3)
                %Estimate snr weights
                snrweights(jx,ix,mx)=abs(psnr(abs(imagesIn{mx}(:,:,jx,ix)),zeros(sizes{1}(1:2))));
                wt(:,:,jx,ix)=snrweights(jx,ix,mx);
            end
        end
    end
    %Calculate mean image
    for kx = [1 2]%Blip up/down images
        imIn=imagesIn{mx}(:,:,:,kx:2:end);
        
        zfsTest(:,:,:,kx) =zfsTest(:,:,:,kx) + sum(wt(:,:,:,kx:2:end).*imIn,4)./abs(max(max(max(sum(imIn,4)))));
    end
    zfsTest=zfsTest./sum(wt(:));
end

%% Calculate and distortion correction
if (peSignConvention) %Even starting frame => +ve peSign
    peSign=[1 -1];
else %Odd starting frame => -ve peSign
    peSign=[-1 1];
end

disp('Calculating distortion correction');
if nargin==4
    save_volume_call_topup(zfsTest, dn, params, peSign);
else
    save_volume_call_topup(zfsTest, dn, params, peSign, varargin{1});
end;
%% Apply to correct other metabolites
disp('Applying correction to other metabolites');
[out{1:nMets, 1}] = deal(zeros(sizes{1}));

if peSignConvention %Even starting frame => +ve peSign
    peSign=mod(1:sizes{1}(4),2)*2-1;
else %Odd starting frame => -ve peSign
    peSign=-mod(1:sizes{1}(4),2)*2+1;
end


for mx = 1:nMets
    out{mx}=apply_computed_topup(imagesIn{mx},dn,params,peSign);
end
