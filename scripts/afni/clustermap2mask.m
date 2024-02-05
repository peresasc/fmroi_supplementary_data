%%
datadir = '/media/andre/data8t/fmroi/fmroi_qc/dataset/afni-clustermask_map';
outdir = '/media/andre/data8t/fmroi/fmroi_qc/dataset/afni-clustermask';
if ~isfolder(outdir)
    mkdir(outdir)
end

roistruc = dir(datadir);
roinames = cell(length(roistruc),1);
for s = 1:length(roistruc)
    if ~roistruc(s).isdir
        roinames{s} = roistruc(s).name;
    end
end
roinames(cellfun(@isempty,roinames)) = [];


for i = 1:length(roinames)
    vmap = spm_vol(fullfile(datadir,roinames{i}));
    roi = spm_data_read(vmap);
    nrois = unique(roi(:));
    nrois(nrois==0) = [];
    
    for n = 1:length(nrois)
        binmask = uint16(roi==nrois(n));

        [~,prefix,~] = fileparts(roinames{i});

        outpath = fullfile(outdir,[prefix,...
            '_cluster_',sprintf('%03d',nrois(n)),'.nii']);
        

        vbin = spm_create_vol(vmap);
        vbin.fname = outpath;
        vbin.pinfo = [1;0;0]; % avoid SPM to rescale the masks
        vbin = spm_write_vol(vbin, binmask);
    end
end