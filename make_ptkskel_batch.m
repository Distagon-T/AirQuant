%% 批量生成 PTKskel：遍历 CT 文件夹与气道掩膜文件夹，按患者名匹配后逐个生成骨架
% =========================================================================
% 逻辑（与之前 python batch 脚本一致）：
%   1. 遍历 CT 输入目录下所有患者子文件夹
%   2. 在每个患者文件夹中找到【层数最多】的 CT .nii.gz（依据 dicom_info.json / niftiinfo）
%   3. 在气道掩膜目录中找到同名患者的 <患者>_airway.nii.gz
%   4. 对每个患者运行 PTKskel 生成（bwskel + 自愈合 + ClinicalAirways 验证）
%   5. 输出到目标目录 <输出>/<患者>_masks/... 并生成 info JSON + 汇总 JSON
% =========================================================================
clear; clc; close all;

% =========================================================================
% 1. 路径配置（按需修改）
% =========================================================================
CT_INPUT_DIR      = 'E:\DICOM\2026-04-nifti';            % 患者 CT 文件夹
AIRWAY_INPUT_DIR  = 'E:\DICOM\2026-04-Airway_out_batch';  % 气道掩膜文件夹
OUTPUT_DIR        = 'E:\DICOM\2026-04-PTKskel';           % 输出目录

% 自愈合内核大小
KERNEL_SIZES = [0, 3, 5, 7];

% =========================================================================
% 2. 收集患者列表（两个目录取交集，按文件夹名匹配）
% =========================================================================
ct_patients = collect_patient_folders(CT_INPUT_DIR);
aw_patients = collect_patient_folders(AIRWAY_INPUT_DIR);
matched = intersect(ct_patients, aw_patients);

fprintf('📦 CT 患者数: %d，气道患者数: %d，匹配成功: %d\n', ...
    numel(ct_patients), numel(aw_patients), numel(matched));

if isempty(matched)
    error('❌ 两个目录没有匹配到任何患者，请检查路径/命名。');
end

mkdir(OUTPUT_DIR);
overall_start = tic;
results = {};
success_count = 0;
skip_count = 0;
fail_count = 0;

% =========================================================================
% 3. 逐个患者处理
% =========================================================================
for idx = 1:numel(matched)
    patient_name = matched{idx};
    fprintf('\n[%d/%d] 正在处理: %s ...\n', idx, numel(matched), patient_name);

    info = struct( ...
        'patient_folder', patient_name, ...
        'ct_dir', fullfile(CT_INPUT_DIR, patient_name), ...
        'airway_dir', fullfile(AIRWAY_INPUT_DIR, patient_name), ...
        'selected_ct', '', ...
        'airway_mask', '', ...
        'skel_output', '', ...
        'status', 'pending', ...
        'elapsed_seconds', [], ...
        'error', [], ...
        'skeleton_voxels', []);

    patient_out_dir = fullfile(OUTPUT_DIR, [patient_name, '_PTKskel']);
    info_json_path = fullfile(patient_out_dir, [patient_name, '_ptkskel_info.json']);
    skel_out = fullfile(patient_out_dir, [patient_name, '_airway_PTKskel.nii.gz']);

    % 【断点续传】
    if exist(fullfile(patient_out_dir, [patient_name, '_airway_PTKskel.nii.gz']), 'file') && ...
       exist(info_json_path, 'file')
        info.status = 'skipped';
        info.skel_output = fullfile(patient_out_dir, [patient_name, '_airway_PTKskel.nii.gz']);
        results{end+1} = info;
        skip_count = skip_count + 1;
        fprintf('   ⏭️ 已存在骨架结果，跳过。\n');
        continue;
    end

    patient_start = tic;
    try
        % --- 3.1 找该患者层数最多的 CT ---
        ct_dir = fullfile(CT_INPUT_DIR, patient_name);
        ct_file = find_largest_slice_nifti(ct_dir);
        assert(~isempty(ct_file), '未找到 CT .nii.gz');
        info.selected_ct = ct_file;
        fprintf('   🎯 选中 CT: %s\n', get_filename(ct_file));

        % --- 3.2 找该患者的气道掩膜 ---
        aw_dir = fullfile(AIRWAY_INPUT_DIR, patient_name);
        aw_file = find_airway_mask(aw_dir, patient_name);
        assert(~isempty(aw_file), '未找到气道掩膜 <患者>_airway.nii.gz');
        info.airway_mask = aw_file;
        fprintf('   🎯 气道掩膜: %s\n', get_filename(aw_file));

        % --- 3.3 加载气道掩膜并清洗 ---
        meta_seg = niftiinfo(aw_file);
        seg_raw  = logical(niftiread(meta_seg));
        seg_base = imfill(seg_raw, 'holes');
        CC = bwconncomp(seg_base, 26);
        numPixels = cellfun(@numel, CC.PixelIdxList);
        [~, idx_max] = max(numPixels);
        seg_clean = false(size(seg_base));
        seg_clean(CC.PixelIdxList{idx_max}) = true;
        seg_base = seg_clean;

        % --- 3.4 自愈合骨架提取 ---
        skel = extract_self_healing_skeleton(seg_base, KERNEL_SIZES, seg_base);
        assert(~isempty(skel), '自愈合骨架提取失败');

        % --- 3.5 保存骨架 ---
        mkdir(patient_out_dir);
        meta_skel = meta_seg;
        meta_skel.Datatype = 'uint8';
        if isfield(meta_skel, 'MultiplicativeScaling'), meta_skel.MultiplicativeScaling = 1; end
        if isfield(meta_skel, 'AdditiveOffset'), meta_skel.AdditiveOffset = 0; end
        niftiwrite(uint8(skel), skel_out, meta_skel, 'Compressed', true);

        info.skel_output = skel_out;
        info.skeleton_voxels = nnz(skel);
        info.elapsed_seconds = toc(patient_start);
        info.status = 'success';

        write_json(info_json_path, info);
        success_count = success_count + 1;
        fprintf('   ✅ 完成！骨架体素数: %d，耗时: %.2f 秒\n', nnz(skel), toc(patient_start));

    catch ME
        info.status = 'failed';
        info.error = ME.message;
        info.elapsed_seconds = toc(patient_start);
        % 失败也写 info json，方便排查
        mkdir(patient_out_dir);
        write_json(info_json_path, info);
        fail_count = fail_count + 1;
        fprintf('   ❌ 失败: %s\n', ME.message);
    end
    results{end+1} = info;
end

% =========================================================================
% 4. 汇总 JSON
% =========================================================================
summary = struct( ...
    'ct_input_dir', CT_INPUT_DIR, ...
    'airway_input_dir', AIRWAY_INPUT_DIR, ...
    'output_dir', OUTPUT_DIR, ...
    'total_patients', numel(matched), ...
    'success', success_count, ...
    'skipped', skip_count, ...
    'failed', fail_count, ...
    'total_elapsed_seconds', toc(overall_start), ...
    'patients', {results});
summary_path = fullfile(OUTPUT_DIR, 'ptkskel_summary.json');
write_json(summary_path, summary);

fprintf('\n======================================================\n');
fprintf('🎉 批量 PTKskel 生成结束！\n');
fprintf('⏱️ 总耗时: %.2f 秒\n', toc(overall_start));
fprintf('📊 统计: 成功 %d，跳过 %d，失败 %d\n', success_count, skip_count, fail_count);
fprintf('💾 汇总: %s\n', summary_path);
fprintf('======================================================\n');

% =========================================================================
% 局部函数
% =========================================================================
function names = collect_patient_folders(base_dir)
    d = dir(base_dir);
    names = {};
    for i = 1:numel(d)
        if d(i).isdir && ~strcmp(d(i).name, '.') && ~strcmp(d(i).name, '..')
            names{end+1} = d(i).name; %#ok<AGROW>
        end
    end
    names = sort(names);
end

function nii = find_largest_slice_nifti(patient_dir)
    % 在患者文件夹中找到层数最多的 .nii.gz（跳过 dicom_info.json）
    files = dir(fullfile(patient_dir, '*.nii.gz'));
    best_file = '';
    best_slices = -1;
    for i = 1:numel(files)
        fpath = fullfile(patient_dir, files(i).name);
        try
            meta = niftiinfo(fpath);
            sz = meta.ImageSize;
            n = sz(3);  % 层数 = z 轴
            if n > best_slices
                best_slices = n;
                best_file = fpath;
            end
        catch
        end
    end
    if best_slices < 0
        nii = '';
    else
        nii = best_file;
    end
end

function aw = find_airway_mask(aw_dir, patient_name)
    % 找到 <患者>_airway.nii.gz（兼容 _airway.nii / .nii.gz）
    aw = '';
    cand = fullfile(aw_dir, [patient_name, '_airway.nii.gz']);
    if exist(cand, 'file')
        aw = cand;
        return;
    end
    % 兜底：目录内任意含 'airway' 的 .nii.gz
    files = dir(fullfile(aw_dir, '*airway*.nii.gz'));
    if ~isempty(files)
        aw = fullfile(aw_dir, files(1).name);
    end
end

function skel = extract_self_healing_skeleton(seg_base, kernel_sizes, ~)
    % 自愈合骨架：kernel 0/3/5/7 逐级净化，直到 ClinicalAirways 能编译
    skel = [];
    for k = kernel_sizes
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
        skel_candidate = bwskel(seg_for_skel, 'MinBranchLength', 15);
        try
            ClinicalAirways(skel_candidate, ...
                'source', double(seg_base), ...
                'seg', seg_base, ...
                'fillholes', 0, 'largestCC', 1, ...
                'plane_sample_sz', 0.5, 'spline_sample_sz', 0.5);
            skel = skel_candidate;
            fprintf('   ✅ 骨架等级 kernel=%d 验证通过\n', k);
            return;
        catch
            fprintf('   ❌ kernel=%d 验证失败\n', k);
        end
    end
end

function write_json(path, data)
    fid = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('无法写入 JSON: %s', path);
    end
    fprintf(fid, '%s', jsonencode(data));
    fclose(fid);
end

function name = get_filename(path)
    [~, name, ext] = fileparts(path);
    name = [name, ext];
end
