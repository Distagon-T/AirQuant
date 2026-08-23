%% 快速测试 skel_2_digraph 的 BFS 断言（不加载大 CT）
function debug_bfs
    function log(s)
        fprintf(2, '[LOG] %s\n', s);
    end

SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';
PAD = 20;

log('S1: load seg');
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));

idx = find(seg_raw);
[i1,i2,i3] = ind2sub(size(seg_raw), idx);
mn = max(1, [min(i1), min(i2), min(i3)] - PAD);
mx = min(size(seg_raw), [max(i1), max(i2), max(i3)] + PAD);
seg_crop = seg_raw(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
clear seg_raw;

seg_base = imfill(seg_crop, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;

log('S2: bwskel');
skel = bwskel(seg_base, 'MinBranchLength', 15);
fprintf(2, '[LOG] skel %d voxels\n', nnz(skel));

log('S3: skel_2_digraph direct');
try
    [g, glink, gnode] = skel_2_digraph(skel, 'topnode');
    fprintf(2, '[LOG] SUCCESS: %d links, %d nodes\n', length(glink), length(gnode));
catch ME
    fprintf(2, '[LOG] FAILED: %s\n', ME.message);
    if strcmp(ME.identifier, 'MATLAB:assertion')
        % 深入检查
        log('S4: inspect graph internals');
        [gadj,gnode2,glink2] = Skel2Graph3D(skel,1);
        G = digraph(gadj);
        bins = conncomp(G);
        fprintf(2, '[LOG] glink=%d nodes=%d bins=%d\n', length(glink2), length(gnode2), max(bins));
    end
end
log('FINISHED');
end
