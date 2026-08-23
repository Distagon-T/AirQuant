%% 最小化崩溃隔离测试
clear; clc; close all;
CT_FILE = 'E:\DICOM\2026-04-nifti\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_4.nii.gz';
SEG_FILE = 'E:\DICOM\2026-04-Airway_out_batch\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323\20130305_CHEN RONG_CT_1.2.840.113619.2.55.3.1678396440.5697.1362438466.323_airway.nii.gz';

disp('STEP1: load CT header');
meta_CT = niftiinfo(CT_FILE);
disp('STEP2: load CT data as double');
source = double(niftiread(meta_CT));
fprintf('   source size: %s, class: %s\n', mat2str(size(source)), class(source));
disp('STEP3: load seg');
meta_seg = niftiinfo(SEG_FILE);
seg_raw = logical(niftiread(meta_seg));
fprintf('   seg size: %s\n', mat2str(size(seg_raw)));
disp('STEP4: clean seg (imfill + largest CC)');
seg_base = imfill(seg_raw, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx_max] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx_max}) = true;
seg_base = seg_clean;
fprintf('   seg voxels: %d\n', nnz(seg_base));
disp('STEP5: bwskel');
skel = bwskel(seg_base, 'MinBranchLength', 15);
fprintf('   skel voxels: %d\n', nnz(skel));
disp('STEP6: ClinicalAirways with source=CT');
AQnet = ClinicalAirways(skel, ...
    'source', source, ...
    'header', meta_CT, ...
    'seg', seg_base, ...
    'fillholes', 0, ...
    'largestCC', 1, ...
    'plane_sample_sz', 0.5, ...
    'spline_sample_sz', 0.5);
disp('STEP7: DONE - AQnet built OK');
fprintf('   tubes: %d, skel size: %s\n', length(AQnet.tubes), mat2str(size(AQnet.skel)));
