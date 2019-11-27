function distortion_corrected=save_volume_call_topup(im, fname, params, peSign, varargin)
% Spit out a HP volume to FSL's topup, let it work its magic, load in the
% resulting nifti image. It is recommended that the input consists of two
% temporally summed volumes, 4D doubles with alternating PE directions in
% the 4th dimension. Likewise, it is suggested that fname is the path to
% the HP fid directory.
%
% After calling topup on high-SNR images, it should be possible to then
% apply the resulting data and distortion field to the bic, lac images.
%
% Inputs: 
% -- im: Input X x Y x Z x PE Sign 4D image (HP) 
% -- fname: directory name of the image. Topup will store files in the fid
% directory. 
% -- params: procpar struct from which things like echo spacing are found 
% -- peSign: vector containing direction of the phase encoding axis 
% Outputs: 
% -- distortion_corrected: loaded in result of topup. 
%
% See also APPLY_COMPUTED_TOPUP
% JM, 2015

if nargin < 4 
    error('Incorrect number of arguments'); 
elseif ndims(im) ~= 4 
    error('Incorrect size of input image'); 
elseif ~isstruct(params) 
    error('Incorrect parameter format. Should be a procpar struct.'); 
elseif ~ischar(fname)
    error('Incorrect filename -- not a string');
elseif length(peSign)~=size(im,4)
    error('Incorrect specification for peSign: length differs to that of the image');
end

% LINCOR_FLAG:
% Set to 1 to do linear algebra; 0 to just apply eddy_correct.
LINCOR_FLAG=1; 

fname=fullfile(fname); 
if ~strcmp(fname(end),filesep) 
    fname(end+1) = filesep; 
end

if size(im,2) ~= size(im,1) %Image has been half FOV'd
 im_out=make_ana(squeeze(im), [params.lro/size(im,1) params.lpe./(2*size(im,2)) params.lpe2./size(im,3)]*10); %Voxel size in mm
else %Trust it
 im_out=make_ana(squeeze(im), [params.lro/size(im,1) params.lpe./(size(im,2)) params.lpe2./size(im,3)]*10); %Voxel size in mm
end
 save_untouch_nii(im_out, [fname '-nifti-original.nii.gz']);

 
fid=fopen([fname '-topup-data.txt'],'w'); 
%fid=fopen(fullfile(fname,'-topup-data.txt'),'w');
if fid<1 
    error('File opening failed');
end

%Note to self: check orientations of 'fsl xy' and 'matlab xy'. 
for i=1:size(im,4)
    %fprintf(fid,'%d %d %d %.4f\n',peSign(i),0, 0, params.esp*params.etl);
    fprintf(fid,'%d %d %d %.4f\n',0,peSign(i), 0, params.esp*params.etl);
end
fclose(fid);

topupcmd=sprintf('cd %s; $FSLDIR/bin/topup -v --imain=%s',fname, ['./-nifti-original']);
topupcmd=sprintf('%s --datain=%s', topupcmd, ['./-topup-data.txt']); 
topupcmd=sprintf('%s --out=%s --fout=%s --iout=%s',topupcmd, ['./-toupout'], ['./-fieldout'], ['./topup-unwarped-images']);



if nargin == 4
    topupcmd=sprintf('%s --subsamp=4,2,1 --fwhm=8,4,0 --miter=15,15,40 --numprec=double --lambda=50,20,2 --ssqlambda=1 --estmov=0 --scale=0',topupcmd);%These work!
    %topupcmd=sprintf('%s --subsamp=4,2,1 --fwhm=8,4,0 --miter=15,15,40 --numprec=double --lambda=50,20,2 --warpres=10,2,2 --ssqlambda=1 --estmov=1 --regmod=membrane_energy --scale=0',topupcmd);
elseif nargin == 5
    if ischar(varargin{1})
        topupcmd=sprintf('%s %s',topupcmd,varargin{1});
        disp('Using custom options string');
        disp(varargin{1});
    else
        error('Check varargin{1} -- should be a string');
    end
end

[status]=call_fsl_local(topupcmd);
if status ~= 0
    error('Calling topup failed'); 
end
disp('Susceptibility corrected; affine registration step'); 
if LINCOR_FLAG==0
eddyCmd=sprintf(' cd %s; rm topup-uneddy-images.ecclog 2> /dev/null;$FSLDIR/bin/eddy_correct ./topup-unwarped-images.nii.gz ./topup-uneddy-images.nii.gz 1;',fname);
else
   eddyCmd=sprintf(' cd %s; rm topup-uneddy-images.ecclog 2> /dev/null;$FSLDIR/bin/eddy_correct ./topup-unwarped-images.nii.gz ./topup-uneddy-images.nii.gz 1;',fname);
   eddyCmd=sprintf('%s mv topup-uneddy-images.ecclog topup-uneddy-1.ecclog; $FSLDIR/bin/eddy_correct ./topup-unwarped-images.nii.gz ./topup-uneddy-images.nii.gz 0;',eddyCmd);
   eddyCmd=sprintf('%s mv topup-uneddy-images.ecclog topup-uneddy-0.ecclog;',eddyCmd);
end
[status]=call_fsl_local(eddyCmd);
if status ~= 0
    error('Calling eddy_correct failed'); 
end

%Assuming topup worked...
% Apply eddy_correct, and register everything to the first value? 

%Note to self: ${FSLDIR}/bin/flirt -in $i -ref ${output}_ref -nosearch -interp ${interpm} -o $i -paddingsize 1 >> ${output}.ecclog
%flirtcmd=sprintf('cd %s; $FSLDIR/bin/flirt -in %s', fname, topup);


distortion_corrected=read_avw([fname 'topup-uneddy-images']); 
 

function [status,output] = call_fsl_local(cmd)
% [status, output] = call_fsl(cmd)
% 
% Wrapper around calls to FSL binaries
% clears LD_LIBRARY_PATH and ensures
% the FSL envrionment variables have been
% set up
% Debian/Ubuntu users should uncomment as
% indicated

%fsldir=getenv('FSLDIR');
fsldir='/usr/local/fsl/';
% Debian/Ubuntu - uncomment the following
%fsllibdir=sprintf('%s/%s', fsldir, 'bin');

if ismac
  dylibpath=getenv('DYLD_LIBRARY_PATH');
  setenv('DYLD_LIBRARY_PATH');
else
  ldlibpath=getenv('LD_LIBRARY_PATH');
  setenv('LD_LIBRARY_PATH');
  % Debian/Ubuntu - uncomment the following
  % setenv('LD_LIBRARY_PATH',fsllibdir);
end

command = sprintf('/bin/bash -c ''. %s/etc/fslconf/fsl.sh; export FSLDIR=%s; export PATH=$PATH:%s/bin/; %s''', fsldir, fsldir, fsldir, cmd);
[status,output] = system(command,'-echo');

if ismac
  setenv('DYLD_LIBRARY_PATH', dylibpath);
else
    setenv('LD_LIBRARY_PATH', ldlibpath);
end

if status
    error('FSL call (%s) failed, %s', command, output)
end
