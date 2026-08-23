%% 精简测试：真实 CT -> ClinicalAirways（-nojvm 下定位崩溃）
function debug_aqnet
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';

log('S1: load CT');
meta_CT = niftiinfo(CT_FILE);
source = double(niftiread(meta_CT));
log(sprintf('S2: source %s', mat2str(size(source))));
log('S3: load seg');
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));
seg_base = imfill(seg_raw, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;
log('S4: bwskel');
skel = bwskel(seg_base, 'MinBranchLength', 15);
log(sprintf('S5: skel %d voxels', nnz(skel)));
log('S6: ClinicalAirways with real CT + header');
AQnet = ClinicalAirways(skel, ...
    'source', source, ...
    'header', meta_CT, ...
    'seg', seg_base, ...
    'fillholes', 0, ...
    'largestCC', 1, ...
    'plane_sample_sz', 0.5, ...
    'spline_sample_sz', 0.5);
log(sprintf('S7: DONE tubes=%d skel=%s', length(AQnet.tubes), mat2str(size(AQnet.skel))));
end
