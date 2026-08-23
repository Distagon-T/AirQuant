%% 对照实验：不传 header 时 ClinicalAirways 是否成功（区分 header 重定向 vs 骨架环）
function debug_nohdr
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
PAD = 20;

log('S1: load');
meta_CT = niftiinfo(CT_FILE);
source_full = double(niftiread(meta_CT));
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));

idx = find(seg_raw);
[i1,i2,i3] = ind2sub(size(seg_raw), idx);
mn = max(1, [min(i1), min(i2), min(i3)] - PAD);
mx = min(size(seg_raw), [max(i1), max(i2), max(i3)] + PAD);
source = source_full(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
seg_crop = seg_raw(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
clear source_full seg_raw;

seg_base = imfill(seg_crop, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;

log('S2: bwskel k0');
skel = bwskel(seg_base, 'MinBranchLength', 15);
fprintf(2, '[LOG] skel %d voxels\n', nnz(skel));

log('S3: ClinicalAirways WITHOUT header, source=seg');
try
    net = ClinicalAirways(skel, ...
        'source', double(seg_base), ...
        'seg', seg_base, ...
        'fillholes', 0, ...
        'largestCC', 1, ...
        'plane_sample_sz', 0.5, ...
        'spline_sample_sz', 0.5);
    fprintf(2, '[LOG] NO-HDR SUCCESS tubes=%d\n', length(net.tubes));
catch ME
    fprintf(2, '[LOG] NO-HDR FAILED: %s\n', ME.message);
end

log('S4: ClinicalAirways WITH header, source=CT');
try
    net2 = ClinicalAirways(skel, ...
        'source', source, ...
        'header', meta_CT, ...
        'seg', seg_base, ...
        'fillholes', 0, ...
        'largestCC', 1, ...
        'plane_sample_sz', 0.5, ...
        'spline_sample_sz', 0.5);
    fprintf(2, '[LOG] HDR SUCCESS tubes=%d\n', length(net2.tubes));
catch ME
    fprintf(2, '[LOG] HDR FAILED: %s\n', ME.message);
end
log('FINISHED');
end
