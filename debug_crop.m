%% 裁剪到 bbox 后再跑 bwskel + ClinicalAirways（验证低内存方案）
function debug_crop
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
PAD = 20;

log('S1: load CT');
meta_CT = niftiinfo(CT_FILE);
source_full = double(niftiread(meta_CT));
log('S2: load seg');
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));

log('S3: compute bbox');
idx = find(seg_raw);
[i1,i2,i3] = ind2sub(size(seg_raw), idx);
mn = [min(i1), min(i2), min(i3)];
mx = [max(i1), max(i2), max(i3)];
sz = size(seg_raw);
mn = max(1, mn - PAD);
mx = min(sz, mx + PAD);
fprintf(2, '[LOG] bbox %s -> %s size %s\n', mat2str(mn), mat2str(mx), mat2str(mx-mn+1));

log('S4: crop source + seg');
source = source_full(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
seg_crop = seg_raw(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
clear source_full seg_raw;
fprintf(2, '[LOG] cropped source %s\n', mat2str(size(source)));

log('S5: clean seg (imfill + largest CC)');
seg_base = imfill(seg_crop, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;
log(sprintf('S6: seg %d voxels', nnz(seg_base)));

log('S7: bwskel');
tic;
skel = bwskel(seg_base, 'MinBranchLength', 15);
fprintf(2, '[LOG] bwskel %.1f s, %d voxels\n', toc, nnz(skel));

log('S8: ClinicalAirways (cropped source + header)');
try
    AQnet = ClinicalAirways(skel, ...
        'source', source, ...
        'header', meta_CT, ...
        'seg', seg_base, ...
        'fillholes', 0, ...
        'largestCC', 1, ...
        'plane_sample_sz', 0.5, ...
        'spline_sample_sz', 0.5);
    fprintf(2, '[LOG] S9: DONE tubes=%d skel=%s lims=%s\n', ...
        length(AQnet.tubes), mat2str(size(AQnet.skel)), mat2str(AQnet.lims));
catch ME
    fprintf(2, '[LOG] S9: ClinicalAirways FAILED: %s\n', ME.message);
end
log('FINISHED');
end
