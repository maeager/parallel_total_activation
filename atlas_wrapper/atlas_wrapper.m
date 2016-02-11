 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts structural atlas into functional space by inverse warping the SPM deformation matrix 
%
% Atlas_wrapper.m uses IBASPM (http://www.thomaskoenig.ch/Lester/ibaspm.htm) and 
% the preprocessing code (Jonas Richiardi, http://miplab.unige.ch/richiardi/software.php)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% segmentFiles{1} = gray matter segmentation spm (c1)
% segmentFiles{2} = white matter segmentation spm (c2)
% segmentFiles{3} = CSF segmentation spm (c3)
% atlasFile = Structural Atlas 'atlas.nii';
% deform = SPM deformation matrix 'seg_sn.mat';
% outputDir = is the path of output atlas to be written... 
% threshold gray matter probability
% new_Auto_Labelling(segmentFiles, atlasFile, deform, outputDir, threshold)
% funFN = functional volume (realigned or mean...)
% atlasStructFN = result of new_Auto_Labelling, structural atlas
% outFN = name of output, atlas in functional space;
% oobList=convert_atlas(funcFN,atlasStructFN,outFN,maxRegionNidx)
