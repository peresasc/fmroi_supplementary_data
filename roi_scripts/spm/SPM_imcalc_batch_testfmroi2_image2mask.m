%==========================================================================
% SPM-imcalc Image Threshold ROI Generation
%
%--------------------------------------------------------------------------
% This script generates Region of Interests (ROIs) based on minimum and
% maximum threshold using the SPM-imcalc tool. It assesses SPM-imcalc's
% functionality and performance in creating ROIs, as reported in the paper
% "fMROI: a simple and adaptable toolbox for easy region-of-interest 
% creation", currently under review.
%
% Author: Daniela Valerio, 2024.
%==========================================================================


%% Min and max have the same value

outdir = '/media/andre/data1/fmROI_tests/MARSBAR/SPM_imcal/syntetic';

[fmroirootdir,~,~] = fileparts(which('fmroi'));
templatesdir = fullfile(fmroirootdir,'etc','test_data');
synpath = fullfile(templatesdir,'complex-shapes.nii.gz');

thresh = [1, 5, 5.2, 5.4, 5.6, 5.8];

for r = 1:length(thresh)
    
    thresh_string = ['i1==',sprintf('%2d',thresh(r))];
    
    
    matlabbatch{1}.spm.util.imcalc.input = {[synpath,',1']};
    matlabbatch{1}.spm.util.imcalc.output = ['spm-img2mask_srcimg_syntetic_threshold_',sprintf('%0.1f',thresh(r)),'_',sprintf('%0.1f',thresh(r)),'.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {'/media/andre/data1/fmROI_tests/MARSBAR/SPM_imcal/syntetic'};
    matlabbatch{1}.spm.util.imcalc.expression = thresh_string;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run',matlabbatch);
    
end


%% Min and max have different values

min = [3 7 9 11];
max= [4 8 10 12];

for r = 1:length(min)
    thresh_string = ['i1>',sprintf('%2d',min(r)),'&','i2<',sprintf('%2d',max(r))];
    
    matlabbatch{1}.spm.util.imcalc.input = {[synpath,',1'];[synpath,',1']};
    matlabbatch{1}.spm.util.imcalc.output = ['spm-img2mask_srcimg_syndata_threshold_',sprintf('%1d',min(r)),'_',sprintf('%1d',max(r)),'.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {'/media/andre/data1/fmROI_tests/MARSBAR/SPM_imcal/syntetic'};
    matlabbatch{1}.spm.util.imcalc.expression = thresh_string;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    
    
    spm_jobman('run',matlabbatch);
    
    
end


%% DMN
dmnpath = fullfile(templatesdir,'default_mode_association-test_z_FDR_0.01.nii.gz');
vsrc = spm_vol(dmnpath);
srcvol = spm_data_read(vsrc);


% min
thr = 16*rand(100,2);
nze = zeros(100,1);


for i = 1:length(thr)
    
    % to have second always bigger than the first
    if thr(i,2) < thr(i,1)
        thr(i,[1,2]) = thr(i,[2,1]);
    end
    
    currdata = srcvol;
    currdata(currdata<thr(i,1)) = 0;
    currdata(currdata>thr(i,2)) = 0;
    
    n = nnz(currdata);
    
    while n == 0
        
        thr(i,:) = 16*rand(1,2);
        if thr(i,2) < thr(i,1)
            thr(i,[1,2]) = thr(i,[2,1]);
        end
        currdata = srcvol;
        currdata(currdata<thr(i,1)) = 0;
        currdata(currdata>thr(i,2)) = 0;
        
        n = nnz(currdata);
    end
    
    nze(i) = n;
    
    thresh_string = ['i1>',sprintf('%2d',thr(i,1)),'&','i2<',sprintf('%2d',thr(i,2))];
    
    
    matlabbatch{1}.spm.util.imcalc.input = {[dmnpath,',1'];[dmnpath,',1']};
                                            
    matlabbatch{1}.spm.util.imcalc.output = ['spm-img2mask_srcimg_dmn_threshold_',sprintf('%1d',thr(i,1)),'_',sprintf('%1d',thr(i,2)),'.nii'];
    
    matlabbatch{1}.spm.util.imcalc.outdir = {'/media/andre/data1/fmROI_tests/MARSBAR/SPM_imcal/dmn'};
    matlabbatch{1}.spm.util.imcalc.expression = thresh_string;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    
    spm_jobman('run',matlabbatch);
    
    
end
