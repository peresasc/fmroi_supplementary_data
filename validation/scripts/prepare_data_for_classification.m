%% load_paths
funcpath = '/media/andre/data8t/ds000030/derivatives/fmroi/funcpath.txt';
funcpath = readcell(funcpath,'Delimiter',[";",'\t']);

confpath = '/media/andre/data8t/ds000030/derivatives/fmroi/confpath.txt';
confpath = readcell(confpath,'Delimiter',[";",'\t']);

participants = readtable('/media/andre/data8t/ds000030/participants.tsv',...
    'FileType','delimitedtext');

%% get_confounds_fmriprep
% Remove derivatives components searching for NaNs. Confounds considereded 
% included white matter signal, global signal, six motion parameters, and 
% both anatomical (aCompCor) and temporal (tCompCor) components.

confsz = zeros(length(confpath),2);
for i = 1:length(confpath)
    curconf = confpath{i};
    [pn,fn,ext] = fileparts(curconf);

    conftable = readtable(curconf,'FileType','delimitedtext');
    conf = table2array(conftable);
    selcols = ~any(isnan(conf));
    conf = conf(:,selcols);
    confsz(i,:) = size(conf);
    writematrix(conf, fullfile(pn,[fn,'_selected.csv']),'Delimiter','tab');
end
disp('DONE!!')

%% check_subject_alignment
% Verifies if each functional image path matches the corresponding confound
% file by comparing subject identifiers embedded in their filenames.
% This step ensures data consistency before further processing.

if length(funcpath) ~= length(confpath)
    error(['The number of functional image paths does not match the number of confounds files.', ...
       newline, 'Please ensure that each subject''s functional data has a corresponding confound file.']);
end

subjcmp = zeros(length(funcpath),1);
for i = 1:length(funcpath)
    curfunc = funcpath{i};
    s1 = strfind(curfunc, '_task-rest_bold');
    a = curfunc(1:s1-1);

    curconf = confpath{i};
    c = strfind(curconf, '_task-rest_bold');
    b = curconf(1:c-1);

    subjcmp(i) = strcmpi(a,b);
end
if any(subjcmp == 0)
    error(['Subject ID mismatch between functional images and confound paths.', ...
        newline, 'Please ensure each functional image is paired with the correct confound based on subject ID.']);
else
    disp('Subject IDs in functional images and confounds are properly matched.');
end

%% extract_participant_info
% Extracts subject IDs from the functional filenames and matches them to 
% the participant information table. Converts diagnosis labels ('CONTROL',
% 'SCHZ','BIPOLAR','ADHD') to numeric codes and stores the resulting
% matrix in a CSV file for use in downstream classification analyses.

funcsubid = zeros(length(funcpath),1);
for i = 1:length(funcpath)
    [~,curfunc,~] = fileparts(funcpath{i});
    s1 = strfind(curfunc,'_task-rest_bold');
    funcsubid(i) = str2double(curfunc(5:s1-1));    
end

% Extract participant IDs and remove the "sub-" prefix
id_strings = participants.participant_id;  % vector of strings, e.g., "sub-00123"
id_numbers = str2double(erase(id_strings, "sub-"));  % convert to numeric IDs

% Map diagnosis labels to numeric codes
diagnosis_map = containers.Map({'CONTROL','SCHZ','BIPOLAR','ADHD'}, 1:4);

% Extract diagnosis strings and convert to numeric codes
diagnosis_strings = participants.diagnosis;  % cell array of char
diagnosis_numbers = cellfun(@(s) diagnosis_map(s),diagnosis_strings);

% Combine IDs and diagnosis codes into an Mx2 matrix
participant_info = [id_numbers, diagnosis_numbers];

[ispres,idx] = ismember(funcsubid, id_numbers);

partinfo = participant_info(idx,:);
writematrix(partinfo,'/media/andre/data8t/ds000030/derivatives/fmroi/participant_info.csv');

%% covert_atlas_to_fmri_shape
% Resamples a high-resolution atlas to match the voxel dimensions and space
% of a target fMRI image using SPM's coregistration tools (optionally with
% FreeSurfer - commented in the end of the script). This step ensures that
% the parcellation atlas is properly aligned with the functional data 
% before ROI extraction or time series analysis.

% Paths
source_nii = '/home/andre/tmp/applymask_test/ds30/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii';       % atlas original 2x2x2mm
target_nii = '/home/andre/tmp/applymask_test/ds30/sub-10159_task-rest_bold_space-MNI152NLin2009cAsym_brainmask.nii';   % imagem fMRI 3x3x4mm
resampled_nii = '/home/andre/tmp/applymask_test/ds30/Schaefer2018_100Parcels_7Networks_3x3x4mm.nii';


% SPM's reslice function
spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.coreg.write.ref = {target_nii};
matlabbatch{1}.spm.spatial.coreg.write.source = {source_nii};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;  % nearest neighbor (preserve labels)
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r_';

spm_jobman('run', matlabbatch);

% The resampled atlas will be saved as 'r_atlas_2mm.nii'
movefile('r_atlas_2mm.nii', resampled_nii);

% Compress to .nii.gz
gzip(resampled_nii);
delete(resampled_nii);

% FreeSurfer transformation
% mri_vol2vol --mov Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz \
% --targ sub-10159_task-rest_bold_space-MNI152NLin2009cAsym_brainmask.nii \
% --o Schaefer2018_100Parcels_7Networks_order_FSLMNI152_fsresamp.nii.gz \
% --regheader --interp nearest