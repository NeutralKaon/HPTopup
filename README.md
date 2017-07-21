# HPTopup
Image-Domain Distortion Correction for Hyperpolarised MRI with an EPI Readout. 

Please see [Susceptibility-induced Distortion Correction in Hyperpolarized Echo Planar Imaging](http://onlinelibrary.wiley.com/doi/10.1002/mrm.26839/pdf) for full instructions about what's going on. 


Introduction 
============
Hyperpolarised experiments are difficult, and those aiming to image cardiac metabolism are subject to a wide variety of distoritons from different sources. When Echo Planar Imaging is used as a rapid imaging readout, susceptibility artefacts usually present themselves as compressions or shearings in the image domain. Many algorithms have been proposed to ameliorate these distortions in the world of proton EPI (e.g. for FMRI), but these often assume that the underlying "truth" is not inherently changing rapidly as a function of time. Moreover, references are not guaranteed to work because the object being imaged may not be readily visible to proton MR.

One alternative approach is to acquire blip-reversed EPI data such that the distortions oscillate between alternate acquisitions, and then use underlying knowledge of the _k_-space trajectory taken to estimate local ∆B0 and hence improve distortions. This strategy also provides information about the degree to which the data were acquired off-resonance, and can hence help co-localise hyperpolarised data and proton references. 

This is the "companion" repository to a paper that explores these issues in full detail, which can also be found [here](Publication-miller-lau-tyler-HPTopup-mrm26839.pdf) If it is useful to you, **please cite it** as: 

```
Miller, JJ; Lau, AZ; Tyler, DJ; Susceptibility-induced distortion correction in hyperpolarized echo planar imaging; Magnetic Resonance in Medicine, July 2017; DOI: 10.1002/mrm.26839. 
```

Installation
============

1. Install fsl, either via http://fsl.fmrib.ox.ac.uk or from your favourite package manager. 
2. Ensure that the FSL libraries are reachable via matlab: you may need to add something like the following to your `startup.m` file: 
```
setenv( 'FSLDIR', '/usr/local/fsl'); %OS-specific
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;
```
3. Download and install the repo, and add everything to your matlab path. 


Example Datasets 
============

The three example datasets discussed in the paper are included (as matlab .MAT files) together with illustrative code.

FAQs. 
============
* Can I use this for proton EPI reconstruction in Matlab? 

Yes, but I have yet to test this. 
