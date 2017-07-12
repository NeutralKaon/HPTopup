function out=apply_computed_topup(im, fname, params, peSign)
%Apply a topup distortion map computed previously, e.g. by
%'save_volume_call_topup'. It is assumed that appropriate topup-produced
%files live in fname/, i.e. fname should be a fid directory.
% JM '15
%
% Inputs:
% im -- 4D matrix of data to correct
% fname -- path/to/file.fid
% params -- varian procpar struct
% peSign -- size(im, 4) vector of phase encoding directions
% Outputs:
% out -- 4D matrix containing supplied images, topup-corrected.
%
% See also SAVE_VOLUME_CALL_TOPUP

if nargin ~= 4
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

%NAUGHTY FLAG:
% Set to 1 to do linear algebra; 0 to just apply eddy_correct.
NAUGHTY_FLAG=1;
%

if ~strcmp(fname(end),filesep)
    fname(end+1)=filesep; 
end

lengthOfIm=size(im,4);
%
im=make_ana(squeeze(im), [params.lro/size(im,1) params.lpe./(2*size(im,2)) params.lpe2./size(im,3)]*10);
save_untouch_nii(im,[fname '-im-applytopup.nii.gz']);

fid=fopen([fname '-applytopup-data.txt'],'w');
if fid<1
    error('File opening failed');
end

%Note to self: check orientations of 'fsl xy' and 'matlab xy'.
for i=1:lengthOfIm;
    %fprintf(fid,'%d %d %d %.4f\n',peSign(i),0, 0, params.esp*params.etl);
    fprintf(fid,'%d %d %d %.4f\n',0,peSign(i), 0, params.esp*params.etl);
end
fclose(fid);


topupcmd=sprintf('cd %s; $FSLDIR/bin/applytopup --verbose --imain=%s',fname, ['-im-applytopup']);
topupcmd=sprintf('%s --datain=%s', topupcmd, ['./-applytopup-data.txt']);
topupcmd=sprintf('%s --topup=%s --out=%s',topupcmd, ['./-toupout'], ['./applytopup-unwarped-images']);
topupcmd=sprintf('%s --inindex=1 --method=jac',topupcmd);


[status]=call_fsl_local(topupcmd);
if status ~= 0
    error('Calling topup failed');
end

%Assuming topup worked...
% Apply eddy_correct, and register everything to the first value?

%Make directory of transformations

mkdir(fname, 'transformations');

%Read in top and bottom.mat, do linear algebra.
if NAUGHTY_FLAG == 1
    constructstring=sprintf('cd %s; head -n 7 ./topup-uneddy-1.ecclog|tail -n 4 > ./transformations/top.mat; tail -n 5 ./topup-uneddy-0.ecclog|head -n 4 > ./transformations/bottom.mat',fname);
    call_fsl_local(constructstring);
    
    top=import_flirt_mat(sprintf('%s/transformations/top.mat',fname));
    bottom=import_flirt_mat(sprintf('%s/transformations/bottom.mat',fname));
    
    roottop=inv(sqrtm(top));
    rootbottom=sqrtm((bottom));
    
   
    outmat=mpower(rootbottom*inv(roottop),0.5);
    %outmat=outmat*roottop;
    %halfway=mpower(m2*inv(m1),0.5)*m1;
    %outmat=(outmat+rootbottom*mpower(inv(rootbottom) * roottop, 0.5))./2; 
    %outmat=m1 * pow( inv( m1 ) * m2, x );
   
    
    
    
    testme=imag(outmat(:))./real(outmat(:));
    if abs(sum(testme(~isnan(testme)))) > 1e-4
        error('Significant imaginary part of transformation matrix: something is likely wrong');
    end
    clear testme;
    
    outmat=real(outmat);
    for i=1:lengthOfIm
        curName=sprintf('%s/transformations/MAT_00%.2d',fname,i-1);
        if mod(i,2)==0
            outmat=sqrtm(bottom); 
            dlmwrite(curName,outmat,'delimiter',' ','precision','%.6f');
        else
            outmat=inv(sqrtm(top));
            dlmwrite(curName,inv(outmat),'delimiter',' ','precision','%.6f');
        end
    end
    
    %Change transformations
    
else
    constructstring=sprintf('cd %s; head -n 7 ./topup-uneddy-images.ecclog|tail -n 4 > ./transformations/top.mat; tail -n 5 ./topup-uneddy-images.ecclog|head -n 4 > ./transformations/bottom.mat',fname);
    call_fsl_local(constructstring);
    
    
    for i=1:lengthOfIm
        constructstring=sprintf('cd %s/transformations/; ',fname);
        if mod(i,2)==1
            constructstring=sprintf('%s cp top.mat MAT_00%.2d',constructstring,i-1);
        else
            constructstring=sprintf('%s cp bottom.mat MAT_00%.2d',constructstring,i-1);
        end
        call_fsl_local(constructstring);
    end
    
end
applyxfm4Dcmd=sprintf('cd %s ; $FSLDIR/bin/applyxfm4D %s ./topup-uneddy-images.nii.gz unwarped-uneddied-output ./transformations/ -fourdigit', fname, './applytopup-unwarped-images');
[status]=call_fsl_local(applyxfm4Dcmd);
if status ~= 0
    error('Calling applyxwfm4D failed');
end


out=read_avw([fname 'unwarped-uneddied-output']);


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
fsldir='/usr/local/fsl';

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

command = sprintf('/bin/sh -c ''. %s/etc/fslconf/fsl.sh; %s''', fsldir, cmd);
[status,output] = system(command,'-echo');

if ismac
    setenv('DYLD_LIBRARY_PATH', dylibpath);
else
    setenv('LD_LIBRARY_PATH', ldlibpath);
end

if status
    error('FSL call (%s) failed, %s', command, output)
end



function top = import_flirt_mat(filename)
% Inport_flirt_mat reads in the /path/to/4x4AffineMatrix.mat returned by
% eddy_correct and topup. It returns the matrix as a matlab matrix.

%% Initialize variables.
delimiter = ' ';

startRow = 1;
endRow = inf;


%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);
%% Create output variable
top = [dataArray{1:end-1}];

