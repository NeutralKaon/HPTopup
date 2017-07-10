# HPTopup
Image-Domain Distortion Correction for Hyperpolarised MRI with an EPI Readout

# Introduction 
Hyperpolarised experiments are difficult, and those aiming to image cardiac metabolism are subject to a wide variety of distoritons from different sources. When Echo Planar Imaging is used as a rapid imaging readout, susceptibility artefacts usually present themselves as compressions or shearings in the image domain. Many algorithms have been proposed to ameliorate these distortions in the world of proton EPI (e.g. for FMRI), but these often assume that the underlying "truth" is not inherently changing rapidly as a function of time. Moreover, references are not guaranteed to work because the object being imaged may not be readily visible to proton MR.

One alternative approach is to acquire blip-reversed EPI data such that the distortions oscillate between alternate acquisitions, and then use underlying knowledge of the _k_-space trajectory taken to estimate local ∆B0 and hence improve distortions. This strategy also provides information about the degree to which the data were acquired off-resonance, and can hence help co-localise hyperpolarised data and proton references. 

This is the "companion" repository to a paper that explores these issues in full detail, and can be found in the "publication" folder. If it is useful to you, **please cite it** as: 

<ComingSoon> 

# Installation

# Example Datasets 

# FAQs. 
