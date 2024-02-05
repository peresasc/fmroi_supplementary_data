%==========================================================================
% SPM-MarsBaR Spherical ROI Generation
%
%--------------------------------------------------------------------------
% This script generates spherical Region of Interests (ROIs) using the
% SPM-MarsBaR software. It assesses SPM-MarsBaR's functionality and
% performance in creating ROIs, as reported in the paper "fMROI: a simple
% and adaptable toolbox for easy region-of-interest creation," currently
% under review.
%
% Author: Daniela Valerio, 2024.
%==========================================================================


% Open marsbar, go to Opition> Edit Options> Base Space for ROIs> Choose T1

root_dir = '/media/andre/data1/fmROI_tests/MARSBAR';
% Directory to store (and load) ROI files
roi_dir = fullfile(root_dir, 'MARSBAR_sphere2');

% MarsBaR version check
if isempty(which('marsbar'))
    error('Need MarsBaR on the path');
end
v = str2num(marsbar('ver'));
if v < 0.35
    error('Batch script only works for MarsBaR >= 0.35');
end
marsbar('on');  % needed to set paths etc

[fmroirootdir,~,~] = fileparts(which('fmroi'));
t1path = fullfile(fmroirootdir,'etc','test_data','T1.nii.gz');

vsrc = spm_vol(t1path);
srcvol = spm_data_read(vsrc);
sig = [-1 1];

point = cell(10,1);

for r = 1:100 % number of rois and size
    
    %Coordinates
    x = (sig(randperm(2,1))*randperm(round(size(srcvol,1)/2-r-3),1))...
        +round(size(srcvol,1)/2);
    y = (sig(randperm(2,1))*randperm(round(size(srcvol,2)/2-r-3),1))...
        +round(size(srcvol,2)/2);
    z = (sig(randperm(2,1))*randperm(round(size(srcvol,3)/2-r-3),1))...
        +round(size(srcvol,3)/2);
    
    point{r} = [x,y,z];
    centre = vsrc.mat*[x,y,z,1]';
%    voxelnatcoordinates = vscr.mat\[x,y,z,1]';
    
    coordinates{r,:} = centre(1:3);
    sphere_centre = coordinates{r}';

    sphere_radius = r;
    area =['marsbar-spheremask_srcimg_t1_radius_' num2str(r) '_center_x' num2str(sphere_centre(1)) 'y' num2str(sphere_centre(2)) 'z' num2str(sphere_centre(3))];
    sphere_roi = maroi_sphere(struct('centre', sphere_centre, 'radius', sphere_radius));
    
   
    % Give it a namefm
    details = [area];
%    details = ['marsbar-spheremask_srcimg_t1_radius_' num2str(r) '_center_x' num2str(voxelnatcoordinates(1)) 'y' num2str(y) 'z' num2str(z)];
    details = details(details ~= ' ');
    roitosave = label(sphere_roi, details);
    
    % save ROI to MarsBaR ROI file, in current directory
    detailsmat = [details, '.mat'];
    saveroi(roitosave, fullfile(roi_dir,detailsmat ));
    
    % also Save as image that can be viewed in MRIcroN
    detailsimg = [details, '.nii'];
    sp = maroi('classdata', 'space'); % default space for images saved as ROIs
    save_as_image(roitosave, fullfile(roi_dir,detailsimg), sp);
    a =['create roi' num2str(r)];
    disp(a)

end
