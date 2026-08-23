%% 批量气道定量：PTKskel 骨架提取 + PI10/内壁/外壁/壁厚/面积测量
% =========================================================================
% 整合 myTest_consolidated.m 的完整流程，按患者批量执行：
%   1. 遍历 CT 目录与气道掩膜目录，按患者文件夹名匹配
%   2. 每个患者：找层数最多的 CT + <患者>_airway.nii.gz
%   3. 自愈合骨架提取 (kernel 0/3/5/7) -> ClinicalAirways 构建 AQnet (source=真实CT)
%   4. 分支标签映射 (AQnet 空间) + FWHM 几何测量 (内径/外径/壁厚/管腔面积/管壁面积)
%   5. 导出 CSV (拓扑+几何融合) + 多标签外壁掩膜 (映射回原始CT空间)
%   6. 计算 Pi10 (线性回归) 并保存回归图 PNG + PDF
%   7. 保存分支可视化图（2D/3D 树状图 + 样条图，可选 3D 表面）PNG + PDF
%   8. 每个患者写 <患者>_airquant_info.json，全部完成写 airquant_summary.json
%   9. 断点续传：CSV+外壁掩膜已存在则跳过
% =========================================================================
clear; clc; close all;

% =========================================================================
% 1. 路径配置（按需修改）
% =========================================================================
CT_INPUT_DIR      = 'E:\DICOM\2026-04-nifti\';             % 患者 CT 文件夹
AIRWAY_INPUT_DIR  = 'E:\DICOM\2026-04-Airway_out_batch\';  % 气道掩膜文件夹
OUTPUT_DIR        = 'E:\DICOM\2026-04-Airway_metrics_tmp\';    % 输出目录

% 测量参数
KERNEL_SIZES   = [0, 3, 5, 7];   % 自愈合断桥内核
WALL_SEARCH_MM = 5;              % 外壁搜索最大距离 (mm)
FWHM_BIN_MM    = 0.2;            % FWHM 距离直方图 bin 宽度
MIN_TUBE_PTS   = 3;              % 少于该骨架点数则忽略该分支
SAVE_PLOTS     = true;           % 是否保存 Pi10 回归图（PNG + PDF）
SAVE_BRANCH_PLOTS = true;       % 是否保存分支可视化图（2D/3D 树状图 + 样条图 + 3D 表面图，PNG + PDF）
SAVE_PLOT3D    = true;          % 是否保存 3D 分割表面 RGB 分节图（资源密集，可按需关闭）
MAX_PATIENTS   = Inf;            % 最多处理多少患者（测试时设为小值，如 2；全量设 Inf）

% =========================================================================
% 2. 收集患者列表（两个目录取交集，按文件夹名匹配）
% =========================================================================
ct_patients = collect_patient_folders(CT_INPUT_DIR);
aw_patients_raw = collect_patient_folders(AIRWAY_INPUT_DIR);
% 气道文件夹名带 '_airway' 后缀，剥离后与 CT 侧匹配
aw_patients = regexprep(aw_patients_raw, '_airway$', '');
matched = intersect(ct_patients, aw_patients);

fprintf('CT 患者数: %d，气道患者数: %d，匹配成功: %d\n', ...
    numel(ct_patients), numel(aw_patients), numel(matched));

if isempty(matched)
    error('两个目录没有匹配到任何患者，请检查路径/命名。');
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
if isfinite(MAX_PATIENTS)
    matched = matched(1:min(MAX_PATIENTS, numel(matched)));
end
for idx = 1:numel(matched)
    patient_name = matched{idx};
    fprintf('\n[%d/%d] 正在处理: %s ...\n', idx, numel(matched), patient_name);

    % 气道文件夹名可能带 '_airway' 后缀（如 2026-07-lung-airway），
    % 也可能不带（如 2026-04-Airway_out_batch），这里自适应判断
    aw_patient_dir = fullfile(AIRWAY_INPUT_DIR, [patient_name, '_airway']);
    if ~exist(aw_patient_dir, 'dir')
        aw_patient_dir = fullfile(AIRWAY_INPUT_DIR, patient_name);
    end

    info = struct( ...
        'patient_folder', patient_name, ...
        'ct_dir', fullfile(CT_INPUT_DIR, patient_name), ...
        'airway_dir', aw_patient_dir, ...
        'selected_ct', '', ...
        'airway_mask', '', ...
        'skel_output', '', ...
        'metrics_csv', '', ...
        'wall_mask', '', ...
        'pi10_png', '', ...
        'pi10_pdf', '', ...
        'tree2d_png', '', ...
        'tree2d_pdf', '', ...
        'tree3d_png', '', ...
        'tree3d_pdf', '', ...
        'spline_png', '', ...
        'spline_pdf', '', ...
        'plot3d_png', '', ...
        'plot3d_pdf', '', ...
        'status', 'pending', ...
        'skeleton_voxels', [], ...
        'num_branches', [], ...
        'num_measured', [], ...
        'Pi10', [], ...
        'elapsed_seconds', [], ...
        'error', []);

    patient_out_dir = fullfile(OUTPUT_DIR, [patient_name, '_airquant']);
    mkdir(patient_out_dir);

    skel_out  = fullfile(patient_out_dir, [patient_name, '_airway_PTKskel.nii.gz']);
    csv_out   = fullfile(patient_out_dir, [patient_name, '_full_metrics.csv']);
    wall_out  = fullfile(patient_out_dir, [patient_name, '_airway_OuterWall.nii.gz']);
    png_out   = fullfile(patient_out_dir, [patient_name, '_pi10.png']);
    pdf_out   = fullfile(patient_out_dir, [patient_name, '_pi10.pdf']);
    info_json_path = fullfile(patient_out_dir, [patient_name, '_airquant_info.json']);

    % 【断点续传】CSV + 外壁掩膜都已存在则跳过
    if exist(csv_out, 'file') && exist(wall_out, 'file') && exist(info_json_path, 'file')
        info.status = 'skipped';
        info.metrics_csv = csv_out;
        info.wall_mask = wall_out;
        info.skel_output = skel_out;
        results{end+1} = info;
        skip_count = skip_count + 1;
        fprintf('   已存在定量结果，跳过。\n');
        continue;
    end

    patient_start = tic;
    try
        % --- 3.1 找该患者层数最多的 CT ---
        ct_file = find_largest_slice_nifti(fullfile(CT_INPUT_DIR, patient_name));
        assert(~isempty(ct_file), '未找到 CT .nii.gz');
        info.selected_ct = ct_file;
        fprintf('   选中 CT: %s\n', get_filename(ct_file));

        % --- 3.2 找该患者的气道掩膜 ---
        aw_file = find_airway_mask(aw_patient_dir, patient_name);
        assert(~isempty(aw_file), '未找到气道掩膜 <患者>_airway.nii.gz');
        info.airway_mask = aw_file;
        fprintf('   气道掩膜: %s\n', get_filename(aw_file));

        % --- 3.3 加载 CT 与掩膜，清洗掩膜 ---
        meta_CT  = niftiinfo(ct_file);
        source   = double(niftiread(meta_CT));
        meta_seg = niftiinfo(aw_file);
        seg_raw  = logical(niftiread(meta_seg));

        spacing = meta_CT.PixelDimensions;
        mean_spacing = mean(spacing);
        voxel_vol = prod(spacing);

        seg_base = imfill(seg_raw, 'holes');
        CC = bwconncomp(seg_base, 26);
        numPixels = cellfun(@numel, CC.PixelIdxList);
        [~, idx_max] = max(numPixels);
        seg_clean = false(size(seg_base));
        seg_clean(CC.PixelIdxList{idx_max}) = true;
        seg_base = seg_clean;

        % --- 3.4 自愈合骨架提取 + 构建 AQnet（source=真实CT） ---
        [skel, AQnet] = build_skel_and_net(seg_base, source, meta_CT, KERNEL_SIZES);
        assert(~isempty(skel), '自愈合骨架提取失败');

        % 保存骨架（与 seg 分开，避免覆盖）
        meta_skel = meta_seg;
        meta_skel.Datatype = 'uint8';
        if isfield(meta_skel, 'MultiplicativeScaling'), meta_skel.MultiplicativeScaling = 1; end
        if isfield(meta_skel, 'AdditiveOffset'), meta_skel.AdditiveOffset = 0; end
        niftiwrite(uint8(skel), skel_out, meta_skel, 'Compressed', true);
        info.skel_output = skel_out;
        info.skeleton_voxels = nnz(skel);
        info.num_branches = length(AQnet.tubes);
        fprintf('   骨架体素数: %d，分支数: %d\n', nnz(skel), length(AQnet.tubes));

        % --- 3.5 分支标签映射（AQnet 空间） ---
        label_matrix = zeros(size(AQnet.skel), 'uint16');
        for i = 1:length(AQnet.tubes)
            label_matrix(AQnet.tubes(i).skelpoints) = i;
        end
        n_labeled = nnz(label_matrix > 0);
        fprintf('   已标记 %d 个骨架点。\n', n_labeled);

        % --- 3.6 FWHM 几何测量 ---
        aq_skel   = AQnet.skel;
        aq_seg    = AQnet.seg;
        aq_source = AQnet.source;
        aq_voxdim = AQnet.voxdim;
        aq_mean_spacing = mean(aq_voxdim);
        aq_voxel_vol = prod(aq_voxdim);

        branches = aq_skel;
        dist_from_lumen = bwdist(aq_seg) * aq_mean_spacing;
        results_geom = zeros(0, 9);
        pad = ceil(WALL_SEARCH_MM / aq_mean_spacing);

        outer_wall_mask = false(size(aq_seg));
        wall_volume_by_branch = zeros(length(AQnet.tubes), 1);

        for k = 1:length(AQnet.tubes)
            idx_pts = find(label_matrix == k);
            if length(idx_pts) < MIN_TUBE_PTS
                continue;
            end

            [d1, d2, d3] = ind2sub(size(aq_seg), idx_pts);
            min1 = max(1, min(d1) - pad); max1 = min(size(aq_seg,1), max(d1) + pad);
            min2 = max(1, min(d2) - pad); max2 = min(size(aq_seg,2), max(d2) + pad);
            min3 = max(1, min(d3) - pad); max3 = min(size(aq_seg,3), max(d3) + pad);

            local_branches = branches(min1:max1, min2:max2, min3:max3);
            local_label    = label_matrix(min1:max1, min2:max2, min3:max3);
            local_seg      = aq_seg(min1:max1, min2:max2, min3:max3);
            local_dist     = dist_from_lumen(min1:max1, min2:max2, min3:max3);
            local_source   = aq_source(min1:max1, min2:max2, min3:max3);

            branch_k_mask = (local_label == k);
            other_branches_mask = (local_branches & (local_label ~= k) & (local_label > 0));

            D_k = bwdist(branch_k_mask);
            if any(other_branches_mask(:))
                D_other = bwdist(other_branches_mask);
                local_M = (D_k < D_other);
            else
                local_M = true(size(local_branches));
            end

            % 内径（由管腔截面积反推）
            L = AQnet.tubes(k).stats.arclength;
            lumen_voxels = sum(local_M(:) & local_seg(:));
            LA = (lumen_voxels * aq_voxel_vol) / L;
            D_in = 2 * sqrt(LA / pi);
            Pi_val = pi * D_in;

            % 管壁 FWHM
            wall_zone = local_M & ~local_seg & (local_dist <= WALL_SEARCH_MM);
            d_vals = local_dist(wall_zone);
            hu_vals = local_source(wall_zone);

            if isempty(d_vals)
                continue;
            end

            edges = 0:FWHM_BIN_MM:WALL_SEARCH_MM;
            bin_centers = edges(1:end-1) + FWHM_BIN_MM/2;
            [~, ~, bin_idx] = histcounts(d_vals, edges);

            mean_hu = zeros(length(bin_centers), 1);
            for b = 1:length(bin_centers)
                if any(bin_idx == b)
                    mean_hu(b) = mean(hu_vals(bin_idx == b));
                else
                    mean_hu(b) = NaN;
                end
            end
            mean_hu = fillmissing(mean_hu, 'linear');

            [peak_hu, peak_idx] = max(mean_hu);
            lung_hu = min(mean_hu(peak_idx:end));
            if isempty(lung_hu) || isnan(lung_hu), lung_hu = -850; end

            half_max = (peak_hu + lung_hu) / 2;
            drop_idx = find(mean_hu(peak_idx:end) < half_max, 1);

            if isempty(drop_idx)
                WT = NaN;
            else
                WT = bin_centers(peak_idx + drop_idx - 1);

                render_WT = max(WT, aq_mean_spacing);
                local_wall_reconstruction = local_M & ~local_seg & (local_dist <= render_WT);
                outer_wall_mask(min1:max1, min2:max2, min3:max3) = ...
                    outer_wall_mask(min1:max1, min2:max2, min3:max3) | local_wall_reconstruction;
                wall_volume_by_branch(k) = nnz(local_wall_reconstruction) * aq_voxel_vol;
            end

            D_out = D_in + 2 * WT;
            WA = pi * (D_out/2)^2 - pi * (D_in/2)^2;
            sqrt_WA = sqrt(WA);
            WA_pct = (WA / (WA + LA)) * 100;

            results_geom = [results_geom; k, LA, WA, WA_pct, D_in, D_out, WT, Pi_val, sqrt_WA];
        end
        fprintf('   完成 %d 根分支的几何测量。\n', size(results_geom, 1));
        info.num_measured = size(results_geom, 1);

        % --- 3.7 导出 CSV（拓扑 + 几何融合） ---
        raw_topo_csv = fullfile(patient_out_dir, [patient_name, '_topology_raw.csv']);
        AQnet.ExportCSV(raw_topo_csv);
        topo_table = readtable(raw_topo_csv);

        VariableNames = {'ID', 'LumenArea_mm2', 'WallArea_mm2', 'WA_pct', ...
                         'Inner_Diameter_mm', 'Outer_Diameter_mm', ...
                         'Wall_Thickness_mm', 'Pi_Perimeter_mm', 'Sqrt_WallArea'};

        if isempty(results_geom)
            geom_table = array2table(zeros(0, 9), 'VariableNames', VariableNames);
        else
            geom_table = array2table(results_geom, 'VariableNames', VariableNames);
        end

        final_table = outerjoin(topo_table, geom_table, 'Keys', 'ID', 'MergeKeys', true);
        writetable(final_table, csv_out);
        info.metrics_csv = csv_out;
        delete(raw_topo_csv);
        fprintf('   CSV 已导出: %s\n', get_filename(csv_out));

        % --- 3.8 导出多标签外壁掩膜（映射回原始 CT 空间） ---
        aff_raw_RAS = meta_CT.Transform.T;
        aff_raw_RAS(4,1:3) = 0;
        aff_raw_RAS = aff_raw_RAS./abs(aff_raw_RAS);
        aff_raw_RAS(isnan(aff_raw_RAS)) = 0;
        aff_raw_LPS = aff_raw_RAS*diag([-1,-1,1,1]);
        aff_LPS_pos_red = abs(aff_raw_LPS(1:3,1:3));
        newaxes = [1,2,3];
        for i = 1:3
            vec = aff_LPS_pos_red(:,i);
            newaxes(i) = find(vec);
        end
        aff_LPS_3x3 = aff_raw_LPS(1:3,1:3);
        aff_LPS_red = zeros(3,3);
        for i = 1:3
            aff_LPS_red(:,i) = aff_LPS_3x3(:,newaxes(i));
        end
        flips = [0,0,0];
        for i = 1:3
            if sum(aff_LPS_red(:,i)) == -1
                flips(i) = 1;
            end
        end

        lims = AQnet.lims;
        sz_re = size(seg_base);
        sz_re = sz_re(newaxes);
        outer_full_re = zeros(sz_re);

        clims = lims;
        for ii = 1:3
            if clims(ii) < 1, clims(ii) = 1; end
        end
        j = 1;
        for ii = 4:6
            if clims(ii) > sz_re(j), clims(ii) = sz_re(j); end
            j = j + 1;
        end

        r1 = clims(1,1):clims(1,2); r2 = clims(2,1):clims(2,2); r3 = clims(3,1):clims(3,2);
        outer_full_re(r1, r2, r3) = outer_wall_mask;

        for i = 1:3
            if flips(i)
                outer_full_re = flip(outer_full_re, i);
            end
        end

        outer_wall_orig = ipermute(outer_full_re, newaxes);
        outer_wall_orig = logical(outer_wall_orig);
        final_multi_label_mask = uint8(seg_base);
        final_multi_label_mask(outer_wall_orig) = 2;

        meta_wall = meta_seg;
        meta_wall.Datatype = 'uint8';
        if isfield(meta_wall, 'MultiplicativeScaling'), meta_wall.MultiplicativeScaling = 1; end
        if isfield(meta_wall, 'AdditiveOffset'), meta_wall.AdditiveOffset = 0; end

        niftiwrite(final_multi_label_mask, wall_out, meta_wall, 'Compressed', true);
        info.wall_mask = wall_out;
        fprintf('   内外壁掩膜已保存: %s (1=内腔, 2=外壁)\n', get_filename(wall_out));

        % --- 3.9 计算 Pi10 + 保存回归图 ---
        valid_idx = ~isnan(final_table.Pi_Perimeter_mm) & ~isnan(final_table.Sqrt_WallArea) & (final_table.Pi_Perimeter_mm > 0);
        X_Pi = final_table.Pi_Perimeter_mm(valid_idx);
        Y_sqrtWA = final_table.Sqrt_WallArea(valid_idx);

        if length(X_Pi) > 5
            p = polyfit(X_Pi, Y_sqrtWA, 1);
            Pi10_val = polyval(p, 10);
            info.Pi10 = Pi10_val;
            fprintf('   Pi10 = %.4f\n', Pi10_val);

            if SAVE_PLOTS
                try
                    fig = figure('Name', 'Pi10 Regression', 'Color', 'w', 'Visible', 'off');
                    scatter(X_Pi, Y_sqrtWA, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'w');
                    hold on;
                    x_fit = linspace(min(X_Pi)*0.8, max(X_Pi)*1.2, 100);
                    plot(x_fit, polyval(p, x_fit), 'r-', 'LineWidth', 2);
                    plot([10 10], [0 Pi10_val], 'k--', 'LineWidth', 1.5);
                    plot([0 10], [Pi10_val Pi10_val], 'k--', 'LineWidth', 1.5);
                    scatter(10, Pi10_val, 80, 'rp', 'filled', 'MarkerEdgeColor', 'k');
                    title(sprintf('Pi10 = %.3f (%s)', Pi10_val, patient_name), 'FontSize', 12, 'FontWeight', 'bold');
                    xlabel('Internal Perimeter (Pi, mm)', 'FontSize', 10, 'FontWeight', 'bold');
                    ylabel('Sqrt Wall Area (mm)', 'FontSize', 10, 'FontWeight', 'bold');
                    grid on;
                    % 保存 PNG（预览用）
                    saveas(fig, png_out);
                    % 保存 PDF（矢量图，适合论文/报告）到同一患者文件夹
                    exportgraphics(fig, pdf_out, 'ContentType', 'vector', 'BackgroundColor', 'white');
                    close(fig);
                    info.pi10_png = png_out;
                    info.pi10_pdf = pdf_out;
                catch
                end
            end
        else
            fprintf('   有效分支过少，无法拟合 Pi10。\n');
        end

        % --- 3.10 保存分支可视化图（AirQuant 官方例子的可视化） ---
        if SAVE_BRANCH_PLOTS
            vis_files = save_branch_visualizations(AQnet, patient_name, ...
                patient_out_dir, SAVE_PLOT3D);
            info.tree2d_png = vis_files.tree2d_png;
            info.tree2d_pdf = vis_files.tree2d_pdf;
            info.tree3d_png = vis_files.tree3d_png;
            info.tree3d_pdf = vis_files.tree3d_pdf;
            info.spline_png = vis_files.spline_png;
            info.spline_pdf = vis_files.spline_pdf;
            info.plot3d_png = vis_files.plot3d_png;
            info.plot3d_pdf = vis_files.plot3d_pdf;
        end

        info.elapsed_seconds = toc(patient_start);
        info.status = 'success';
        write_json(info_json_path, info);
        success_count = success_count + 1;
        fprintf('   完成！耗时 %.2f 秒\n', toc(patient_start));

    catch ME
        info.status = 'failed';
        info.error = ME.message;
        info.elapsed_seconds = toc(patient_start);
        mkdir(patient_out_dir);
        write_json(info_json_path, info);
        fail_count = fail_count + 1;
        fprintf('   失败: %s\n', ME.message);
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
summary_path = fullfile(OUTPUT_DIR, 'airquant_summary.json');
write_json(summary_path, summary);

fprintf('\n======================================================\n');
fprintf('批量气道定量结束！\n');
fprintf('总耗时: %.2f 秒\n', toc(overall_start));
fprintf('统计: 成功 %d，跳过 %d，失败 %d\n', success_count, skip_count, fail_count);
fprintf('汇总: %s\n', summary_path);
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
    files = dir(fullfile(patient_dir, '*.nii.gz'));
    best_file = '';
    best_slices = -1;
    for i = 1:numel(files)
        fpath = fullfile(patient_dir, files(i).name);
        try
            meta = niftiinfo(fpath);
            n = meta.ImageSize(3);
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
    aw = '';
    cand = fullfile(aw_dir, [patient_name, '_airway.nii.gz']);
    if exist(cand, 'file')
        aw = cand;
        return;
    end
    files = dir(fullfile(aw_dir, '*airway*.nii.gz'));
    if ~isempty(files)
        aw = fullfile(aw_dir, files(1).name);
    end
end

function [skel, AQnet] = build_skel_and_net(seg_base, source, meta_CT, kernel_sizes)
    skel = [];
    AQnet = [];
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
            net = ClinicalAirways(skel_candidate, ...
                'source', source, ...
                'header', meta_CT, ...
                'seg', seg_base, ...
                'fillholes', 0, ...
                'largestCC', 1, ...
                'plane_sample_sz', 0.5, ...
                'spline_sample_sz', 0.5);
            skel = skel_candidate;
            AQnet = net;
            fprintf('   骨架等级 kernel=%d 验证通过\n', k);
            return;
        catch ME
            fprintf('   kernel=%d 验证失败: %s\n', k, ME.message);
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

function files = save_branch_visualizations(AQnet, patient_name, patient_out_dir, want_plot3d)
    % 保存 AirQuant 官方例子的分支可视化图（PNG + PDF 两种格式）
    % 输出: 2D 树状图 / 3D 树状图 / 3D 样条图 / 3D 分割表面图（均按 RGB 彩色分节着色）
    files = struct('tree2d_png','', 'tree2d_pdf','', ...
                   'tree3d_png','', 'tree3d_pdf','', ...
                   'spline_png','', 'spline_pdf','', ...
                   'plot3d_png','', 'plot3d_pdf','');

    % 尝试肺叶分类（用于 lobe 彩色着色）；失败则改用代数(generation)着色，
    % 保证始终输出 RGB 彩色分节图，而不是黑色单色图
    try
        AQnet.ClassifyLungLobes();
        colour_field = 'lobe';
    catch
        colour_field = 'generation';
        fprintf('   肺叶分类失败，改用代数(generation) RGB 着色。\n');
    end

    % --- 2D 树状图（按 lobe/generation 彩色分节） ---
    try
        fig = figure('Name', 'Airway Tree 2D', 'Color', 'w', 'Visible', 'off');
        AQnet.Plot('colour', colour_field, 'weight', 'generation', 'weightfactor', 10);
        title(sprintf('Airway Tree 2D (%s)', patient_name), 'FontSize', 12, 'FontWeight', 'bold');
        files.tree2d_png = fullfile(patient_out_dir, [patient_name, '_tree2d.png']);
        files.tree2d_pdf = fullfile(patient_out_dir, [patient_name, '_tree2d.pdf']);
        exportgraphics(gcf, files.tree2d_png, 'Resolution', 200, 'BackgroundColor', 'white');
        exportgraphics(gcf, files.tree2d_pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
        close(fig);
        fprintf('   2D 树状图已保存\n');
    catch ME
        fprintf('   2D 树状图保存失败: %s\n', ME.message);
    end

    % --- 3D 树状图（按 lobe/generation 彩色分节） ---
    try
        fig = figure('Name', 'Airway Tree 3D', 'Color', 'w', 'Visible', 'off');
        AQnet.Plot3('colour', colour_field);
        title(sprintf('Airway Tree 3D (%s)', patient_name), 'FontSize', 12, 'FontWeight', 'bold');
        files.tree3d_png = fullfile(patient_out_dir, [patient_name, '_tree3d.png']);
        files.tree3d_pdf = fullfile(patient_out_dir, [patient_name, '_tree3d.pdf']);
        exportgraphics(gcf, files.tree3d_png, 'Resolution', 200, 'BackgroundColor', 'white');
        exportgraphics(gcf, files.tree3d_pdf, 'BackgroundColor', 'white');
        close(fig);
        fprintf('   3D 树状图已保存\n');
    catch ME
        fprintf('   3D 树状图保存失败: %s\n', ME.message);
    end

    % --- 3D 样条图（按 lobe/generation 彩色分节） ---
    try
        fig = figure('Name', 'Airway Splines', 'Color', 'w', 'Visible', 'off');
        AQnet.PlotSpline('colour', colour_field);
        title(sprintf('Airway Splines (%s)', patient_name), 'FontSize', 12, 'FontWeight', 'bold');
        files.spline_png = fullfile(patient_out_dir, [patient_name, '_spline.png']);
        files.spline_pdf = fullfile(patient_out_dir, [patient_name, '_spline.pdf']);
        exportgraphics(gcf, files.spline_png, 'Resolution', 200, 'BackgroundColor', 'white');
        exportgraphics(gcf, files.spline_pdf, 'BackgroundColor', 'white');
        close(fig);
        fprintf('   3D 样条图已保存\n');
    catch ME
        fprintf('   3D 样条图保存失败: %s\n', ME.message);
    end

    % --- 3D 分割表面图（按 lobe/generation RGB 彩色分节） ---
    if want_plot3d
        try
            fig = figure('Name', 'Airway 3D Surface', 'Color', 'w', 'Visible', 'off');
            AQnet.Plot3D('colour', colour_field);
            title(sprintf('Airway 3D Surface (%s)', patient_name), 'FontSize', 12, 'FontWeight', 'bold');
            files.plot3d_png = fullfile(patient_out_dir, [patient_name, '_plot3d.png']);
            files.plot3d_pdf = fullfile(patient_out_dir, [patient_name, '_plot3d.pdf']);
            exportgraphics(gcf, files.plot3d_png, 'Resolution', 200, 'BackgroundColor', 'white');
            exportgraphics(gcf, files.plot3d_pdf, 'BackgroundColor', 'white');
            close(fig);
            fprintf('   3D 分割表面图已保存\n');
        catch ME
            fprintf('   3D 分割表面图保存失败: %s\n', ME.message);
        end
    end
end
