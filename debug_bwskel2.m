%% 用 stderr 逐步定位崩溃（MATLAB 崩溃时 stderr 也会保留）
function debug_bwskel2
    function log(s)
        fprintf(2, '[LOG] %s\n', s);   % 写 stderr，崩溃时也保留
    end

SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
log('STEP1: load seg');
meta_seg = niftiinfo(SEG_FILE);
log('STEP1b: niftiread');
seg_raw = logical(niftiread(meta_seg));
log('STEP2: imfill');
seg_base = imfill(seg_raw, 'holes');
log('STEP3: bwconncomp 26');
CC = bwconncomp(seg_base, 26);
log('STEP4: largest CC');
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;
log(sprintf('seg voxels: %d', nnz(seg_base)));
log('STEP5: bwskel');
tic;
skel = bwskel(seg_base, 'MinBranchLength', 15);
log(sprintf('bwskel done %.2f s, voxels %d', toc, nnz(skel)));
log('DONE');
end
