%% 完整复刻 build_skel_and_net 的 kernel 循环，定位硬崩溃点
function debug_aqnet2
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
KERNEL_SIZES = [0, 3, 5, 7];

log('S1: load CT');
meta_CT = niftiinfo(CT_FILE);
source = double(niftiread(meta_CT));
log('S2: load seg');
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));
log('S3: clean seg');
seg_base = imfill(seg_raw, 'holes');
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
    log(sprintf('K%d: preparing seg_for_skel', k));
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
    log(sprintf('K%d: skel %d voxels, calling ClinicalAirways', k, nnz(skel_candidate)));
    try
        net = ClinicalAirways(skel_candidate, ...
            'source', source, ...
            'header', meta_CT, ...
            'seg', seg_base, ...
            'fillholes', 0, ...
            'largestCC', 1, ...
            'plane_sample_sz', 0.5, ...
            'spline_sample_sz', 0.5);
        log(sprintf('K%d: SUCCESS tubes=%d', k, length(net.tubes)));
        skel = skel_candidate;
        AQnet = net;
        break;
    catch ME
        log(sprintf('K%d: FAILED: %s', k, ME.message));
    end
end
log(sprintf('DONE: skel empty=%d', isempty(skel)));
end
