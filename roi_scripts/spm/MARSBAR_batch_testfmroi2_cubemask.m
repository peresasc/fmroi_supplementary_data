%==========================================================================
% SPM-MarsBaR Cubic ROI Generation
%
%--------------------------------------------------------------------------
% This script generates cubic Region of Interests (ROIs) using the
% SPM-MarsBaR software. It assesses SPM-MarsBaR's functionality and
% performance in creating ROIs, as reported in the paper "fMROI: a simple
% and adaptable toolbox for easy region-of-interest creation," currently
% under review.
%
% Author: Daniela Valerio, 2024.
%==========================================================================

% Open marsbar, go to Opition> Edit Options> Base Space for ROIs> Choose T1
% OR open the varibale MARS.mat no workspace

% Start marsbar to make sure spm_get works
marsbar('on')

root_dir = '/media/andre/data1/fmROI_tests/MARSBAR';
% Directory to store (and load) ROI files
roi_dir = fullfile(root_dir, 'MARSBAR_cube3');

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
  
    x = (sig(randperm(2,1))*randperm(round(size(srcvol,1)/2-r/2-3),1))...
        +round(size(srcvol,1)/2);
    y = (sig(randperm(2,1))*randperm(round(size(srcvol,2)/2-r/2-3),1))...
        +round(size(srcvol,2)/2);
    z = (sig(randperm(2,1))*randperm(round(size(srcvol,3)/2-r/2-3),1))...
        +round(size(srcvol,3)/2);
    
    point{r} = [x,y,z];
    centre = vsrc.mat*[x,y,z,1]';
    
    coordinates{r,:} = centre(1:3);
    cube_centre = coordinates{r}';
    
   
    if (-1)^r == 1 % even
        cube_centre_positive = cube_centre + r/2;
        cube_centre_negative = cube_centre - r/2;
        cube_centre_total = [cube_centre_negative; cube_centre_positive];
        box_widths = abs(diff(cube_centre_total));
        
    else % odd
        
        cube_centre_positive = cube_centre + fix(r/2);
        cube_centre_negative = cube_centre - fix(r/2);
        cube_centre_total = [cube_centre_negative; cube_centre_positive];
        box_widths = abs(diff(cube_centre_total))+1;
        
    end
    

    area =['spm-cubicmask_srcimg_t1_edge_',sprintf('%03d',r),'_center_x',...
        num2str(cube_centre(1)),'y',num2str(cube_centre(2)),'z',num2str(cube_centre(3))];
    cube_roi = maroi_box(struct('centre', cube_centre, 'widths', box_widths));
    
   
    % Give it a name
    details = ['spm-cubicmask_srcimg_t1_edge_',sprintf('%03d',r),'_center_x',...
        num2str(x),'y',num2str(y),'z',num2str(z)];
    details = details(details ~= ' ');
    roitosave = label(cube_roi, details);
    
    % save ROI to MarsBaR ROI file, in current directory
    detailsmat = [details, '.mat'];
    if ~isfolder(roi_dir)
        mkdir(roi_dir)
    end
    saveroi(roitosave, fullfile(roi_dir,detailsmat ));
    
    % also Save as image that can be viewed in MRIcroN
    detailsimg = [details, '.nii'];
    sp = maroi('classdata', MARS.OPTIONS.spacebase.fname);
    save_as_image(roitosave, fullfile(roi_dir,detailsimg), sp);
    a =['create roi' num2str(r)];
    disp(a)

end