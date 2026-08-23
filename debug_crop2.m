%% 裁剪 + kernel 自愈合循环验证（完整组合）
function debug_crop2
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
KERNEL_SIZES = [0, 3, 5, 7];
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
mn = max(1, [min(i1), min(i2), min(i3)] - PAD);
mx = min(size(seg_raw), [max(i1), max(i2), max(i3)] + PAD);
fprintf(2, '[LOG] bbox %s -> %s\n', mat2str(mn), mat2str(mx));

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
log(sprintf('S4: seg %d voxels', nnz(seg_base)));

skel = [];
AQnet = [];
for k = KERNEL_SIZES
    log(sprintf('K%d: prep', k));
    if k == 0
        seg_for_skel = seg_base;
    else
        se = strel('cube', k);
        seg_for_skel = imopen(seg_base, se);
        CC_skel = bwconncomp(seg_for_skel, 26);
        numPixels_skel = cellfun(@numel, CC_skel.PixelIdxList);
        [~, idx_skel] = max(numPixels_skel);
        temp = false(size(seg_for_skel));
        temp(CC_skel.PixelIdxList{idx_skel}) = true;
        seg_for_skel = temp;
    end
    log(sprintf('K%d: bwskel', k));
    skel_candidate = bwskel(seg_for_skel, 'MinBranchLength', 15);
    fprintf(2, '[LOG] K%d skel %d voxels, ClinicalAirways...\n', k, nnz(skel_candidate));
    try
        net = ClinicalAirways(skel_candidate, ...
            'source', source, ...
            'header', meta_CT, ...
            'seg', seg_base, ...
            'fillholes', 0, ...
            'largestCC', 1, ...
            'plane_sample_sz', 0.5, ...
            'spline_sample_sz', 0.5);
        fprintf(2, '[LOG] K%d SUCCESS tubes=%d lims=%s\n', k, length(net.tubes), mat2str(net.lims));
        skel = skel_candidate;
        AQnet = net;
        break;
    catch ME
        fprintf(2, '[LOG] K%d FAILED: %s\n', k, ME.message);
    end
end
log(sprintf('DONE: skel empty=%d', isempty(skel)));
end
