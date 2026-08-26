%% 气道重塑定量特征扩展：形态学分级 + 拓扑网络 + 管壁密度纹理/边界模糊 + T/D 急剧变化 + PCA
% =========================================================================
% 读取 batch_airway_quant.m 生成的 _full_metrics.csv（AirQuant 拓扑+几何），
% 并重新加载 CT/掩膜做 FWHM 提取管壁 HU 密度/纹理与边界模糊度统计。
% 输出：每患者一行的特征表 airway_features.csv + 汇总（便于后续建模/PCA）
%
% 四大类特征（对齐“跨病种、泛气道急性炎性重塑 Pan-airway Acute Inflammatory Remodeling”）：
%   一、形态学分级：各代 WA% 均值、T/D (WT/OD)、内径/外径、支气管扩张参考(BAR 需动脉)
%       + T/D 急剧变化：跨树 std/CV、沿代数稳健斜率、近端→远端阶跃、局灶 |z|>2 异常分支占比
%   二、拓扑网络：分叉角(parent/sibling)、终端分支数/修剪率、迂曲度、分支数、最大代数
%   三、管壁密度纹理：Wall Mean HU / std / skew / kurt（全树+分代）+ PCA 异质性
%   四、FWHM 边界模糊度：壁峰值 HU / 管壁-肺对比度 / 外缘过渡带宽度 / 边界锐度(HU/mm)
%       （炎症水肿 → HU 下降、对比度降低、过渡带变宽、锐度减小），并按 >=5 代聚合
% =========================================================================
clear; clc; close all;

% =========================================================================
% 1. 路径配置
% =========================================================================
% 路径通过环境变量配置（避免硬编码个人路径）；未设置时使用当前目录下的默认子目录
METRICS_DIR   = getenv('AQ_METRICS_DIR');   % batch_airway_quant 输出目录
if isempty(METRICS_DIR)
    METRICS_DIR = fullfile(pwd, 'metrics');
end
FEATURES_DIR  = getenv('AQ_FEATURES_DIR');  % 本脚本输出目录
if isempty(FEATURES_DIR)
    FEATURES_DIR = fullfile(pwd, 'features');
end

% 测量参数（与 batch_airway_quant 保持一致）
KERNEL_SIZES   = [0, 3, 5, 7];
WALL_SEARCH_MM = 5;
FWHM_BIN_MM    = 0.2;
MIN_TUBE_PTS   = 3;
% 最多处理多少患者：环境变量 AQ_MAX_PATIENTS 可覆盖（小样本试跑用，全量设 Inf）
MAX_PATIENTS = str2double(getenv('AQ_MAX_PATIENTS'));
if isnan(MAX_PATIENTS) || MAX_PATIENTS <= 0
    MAX_PATIENTS = Inf;
end
% 是否计算第三类【管壁密度/纹理 + PCA】。
% 需要重读 CT + 掩膜做 FWHM，每个患者约 1-2 分钟。
% 已在 MATLAB GUI 里跑（内存充足），建议保持 true 以输出全部三大类特征。
COMPUTE_DENSITOMETRY = true;

mkdir(FEATURES_DIR);

% =========================================================================
% 2. 收集患者（METRICS_DIR 下的 *_airquant 目录）
% =========================================================================
subdirs = dir(fullfile(METRICS_DIR, '*_airquant'));
patient_names = {};
for i = 1:numel(subdirs)
    if subdirs(i).isdir && subdirs(i).name(end) ~= '.'
        patient_names{end+1} = subdirs(i).name(1:end-9);   % 去掉末尾 '_airquant'
    end
end
patient_names = sort(patient_names);
if isfinite(MAX_PATIENTS)
    patient_names = patient_names(1:min(MAX_PATIENTS, numel(patient_names)));
    fprintf('小样本试跑模式：仅处理 %d 个患者（AQ_MAX_PATIENTS=%d）。\n', ...
        numel(patient_names), MAX_PATIENTS);
end
fprintf('找到 %d 个患者定量结果。\n', numel(patient_names));

% =========================================================================
% 3. 逐个患者计算特征
% =========================================================================
feature_table = table;
for idx = 1:numel(patient_names)
    patient = patient_names{idx};
    fprintf('\n[%d/%d] 计算特征: %s ...\n', idx, numel(patient_names), patient);
    try
    pdir = fullfile(METRICS_DIR, [patient, '_airquant']);
    csv_file = fullfile(pdir, [patient, '_full_metrics.csv']);
    info_file = fullfile(pdir, [patient, '_airquant_info.json']);

    % 读取 info JSON（拿 CT / 掩膜路径 和 Pi10）
    ct_file = ''; seg_file = ''; skel_file = ''; Pi10 = NaN;
    if exist(info_file, 'file')
        info = jsondecode(fileread(info_file));
        if isfield(info, 'selected_ct'),  ct_file  = info.selected_ct;  end
        if isfield(info, 'airway_mask'),  seg_file = info.airway_mask;  end
        % AirQuant 已生成的骨架（batch_airway_quant 保存，拓扑已验证）
        % densitometry 优先用它构建 ClinicalAirways，避免 bwskel 生成异常拓扑
        if isfield(info, 'skel_output') && ~isempty(info.skel_output)
            skel_file = info.skel_output;
        end
        % Pi10 强制为标量 double（JSON 中可能为 []/null/字符串，否则表拼接类型不一致崩溃）
        if isfield(info, 'Pi10')
            p10 = info.Pi10;
            if ischar(p10) || isstring(p10)
                p10 = str2double(p10);
            end
            if isnumeric(p10) && ~isempty(p10)
                Pi10 = p10(1);
            end
        end
    end

    % --- 读取 CSV（AirQuant 拓扑 + 几何） ---
    if ~exist(csv_file, 'file')
        fprintf('   跳过（无 CSV）: %s\n', csv_file);
        continue;
    end
    T = readtable(csv_file);
    if isempty(T) || height(T) == 0
        fprintf('   跳过（空 CSV）\n');
        continue;
    end

    % 确保必要列存在
    need = {'generation','WA_pct','Wall_Thickness_mm','Outer_Diameter_mm', ...
            'Inner_Diameter_mm','LumenArea_mm2','WallArea_mm2','Pi_Perimeter_mm', ...
            'stats_tortuosity','stats_parent_deg','stats_sibling_deg'};
    for n = need
        if ~ismember(n{1}, T.Properties.VariableNames)
            T.(n{1}) = nan(height(T), 1);
        end
    end

    gen = T.generation;
    WA  = T.WA_pct;
    WT  = T.Wall_Thickness_mm;
    OD  = T.Outer_Diameter_mm;
    Din = T.Inner_Diameter_mm;
    LA  = T.LumenArea_mm2;
    WAA = T.WallArea_mm2;
    Pi  = T.Pi_Perimeter_mm;
    tort = T.stats_tortuosity;
    pdeg = T.stats_parent_deg;
    sdeg = T.stats_sibling_deg;
    % children 列（判断终端分支）
    hasch = zeros(height(T),1);
    if ismember('children_1', T.Properties.VariableNames)
        hasch = ~isnan(T.children_1);
    end
    if ismember('children_2', T.Properties.VariableNames)
        hasch = hasch | ~isnan(T.children_2);
    end

    % =====================================================================
    % 特征占位
    % =====================================================================
    f = struct();
    f.patient_folder = patient;
    f.Pi10 = Pi10;

    % ---- 一、形态学分级特征 ----
    f.n_branches      = height(T);
    f.Din_mean_all    = nanmean(Din);
    f.Din_mean_gen3   = nanmean(Din(gen==3));
    f.Din_mean_gen4   = nanmean(Din(gen==4));
    f.Din_mean_gen5   = nanmean(Din(gen==5));
    f.Dout_mean_all   = nanmean(OD);
    % 各代 WA% 均值（远端小气道敏感）
    f.WA_pct_gen3     = nanmean(WA(gen==3));
    f.WA_pct_gen4     = nanmean(WA(gen==4));
    f.WA_pct_gen5     = nanmean(WA(gen==5));
    f.WA_pct_gen3to6  = nanmean(WA(gen>=3 & gen<=6));
    f.WA_pct_all      = nanmean(WA);
    % T/D = WT / OD（无量纲）
    TDr = WT ./ OD;
    TDr(~isfinite(TDr) | OD <= 0) = NaN;   % 过滤 OD<=0 / Inf 异常分支
    f.TD_ratio_all    = nanmean(TDr);
    f.TD_ratio_gen3   = nanmean(TDr(gen==3));
    f.TD_ratio_gen4   = nanmean(TDr(gen==4));
    f.TD_ratio_gen5   = nanmean(TDr(gen==5));
    % ---- T/D 急剧变化（急性炎性重塑：跨树异质性 + 沿代数陡变 + 局灶异常）----
    f.TD_ratio_std_all      = nanstd(TDr);
    if f.TD_ratio_all > 0
        f.TD_ratio_cv_all = f.TD_ratio_std_all / f.TD_ratio_all;   % 变异系数
    else
        f.TD_ratio_cv_all = NaN;
    end
    f.TD_ratio_std_gen5plus = nanstd(TDr(gen>=5));
    f.TD_slope_vs_gen       = td_slope_vs_gen(gen, TDr);          % T/D 随代数稳健斜率（每代变化量）
    f.TD_outlier_ratio_z2   = td_outlier_ratio(TDr);              % |z|>2 局灶陡变分支占比
    f.TD_distal_minus_proximal = nanmean(TDr(gen>=5)) - nanmean(TDr(gen<=4)); % 近端→远端阶跃
    % 管腔/管壁面积
    f.LA_mean_all     = nanmean(LA);
    f.WA_mean_all     = nanmean(WAA);
    f.Pi_mean_all     = nanmean(Pi);
    % 支气管扩张参考：内径均值（BAR 需伴行动脉掩膜，暂输出 Din 供后续计算）

    % ---- 二、拓扑网络特征 ----
    f.max_generation      = max(gen);
    f.mean_tortuosity     = nanmean(tort);
    f.std_tortuosity      = nanstd(tort);
    f.mean_parent_angle   = nanmean(pdeg);     % 母-子分叉角
    f.mean_sibling_angle  = nanmean(sdeg);     % 兄弟分支角
    f.mean_parent_angle_gen3 = nanmean(pdeg(gen==3));
    f.mean_parent_angle_gen4 = nanmean(pdeg(gen==4));
    % 终端分支（无子代）数量与修剪率
    f.n_terminal_total    = sum(~hasch);
    f.n_terminal_gen5plus = sum((gen>=5) & ~hasch);
    f.n_terminal_gen6plus = sum((gen>=6) & ~hasch);
    f.pruning_ratio_gen5  = f.n_terminal_gen5plus / max(1, f.n_branches);   % 高代终端占比
    f.pruning_ratio_gen6  = f.n_terminal_gen6plus / max(1, f.n_branches);
    f.mean_WA_pct_terminal = nanmean(WA(~hasch));

    % ---- 三、管壁密度/纹理（需重读 CT + FWHM） ----
    f.wall_hu_mean = NaN; f.wall_hu_std = NaN; f.wall_hu_skew = NaN; f.wall_hu_kurt = NaN;
    f.wall_hu_mean_gen3 = NaN; f.wall_hu_mean_gen4 = NaN; f.wall_hu_mean_gen5 = NaN;
    f.pca_explained_1 = NaN; f.pca_explained_2 = NaN; f.pca_explained_3 = NaN;
    f.pca_first_pc_std = NaN;
    % ---- 四、FWHM 边界模糊度 + FWHM 版 T/D 占位（由 compute_wall_densitometry 填充）----
    f.blur_peak_hu_mean = NaN; f.blur_peak_hu_std = NaN; f.blur_lung_hu_mean = NaN;
    f.blur_contrast_mean = NaN; f.blur_contrast_std = NaN;
    f.blur_trans_width_mean = NaN; f.blur_trans_width_std = NaN;
    f.blur_edge_sharp_mean = NaN; f.blur_edge_sharp_std = NaN;
    f.blur_contrast_gen5plus = NaN; f.blur_trans_width_gen5plus = NaN;
    f.blur_edge_sharp_gen5plus = NaN; f.blur_peak_hu_gen5plus = NaN;
    f.TD_fwhm_all = NaN; f.TD_fwhm_std = NaN; f.TD_fwhm_cv = NaN;
    f.TD_fwhm_gen5plus = NaN; f.TD_fwhm_std_gen5plus = NaN; f.TD_fwhm_slope_vs_gen = NaN;

    if COMPUTE_DENSITOMETRY && ~isempty(ct_file) && ~isempty(seg_file)
        try
            wall_stats = compute_wall_densitometry(ct_file, seg_file, skel_file, ...
                KERNEL_SIZES, WALL_SEARCH_MM, FWHM_BIN_MM, MIN_TUBE_PTS);
            if ~isempty(wall_stats)
                f.wall_hu_mean = wall_stats.hu_mean_all;
                f.wall_hu_std  = wall_stats.hu_std_all;
                f.wall_hu_skew = wall_stats.hu_skew_all;
                f.wall_hu_kurt = wall_stats.hu_kurt_all;
                f.wall_hu_mean_gen3 = wall_stats.hu_mean_gen3;
                f.wall_hu_mean_gen4 = wall_stats.hu_mean_gen4;
                f.wall_hu_mean_gen5 = wall_stats.hu_mean_gen5;
                f.pca_explained_1 = wall_stats.pca_explained(1);
                f.pca_explained_2 = wall_stats.pca_explained(2);
                f.pca_explained_3 = wall_stats.pca_explained(3);
                f.pca_first_pc_std = wall_stats.pca_pc1_std;
                % FWHM 边界模糊度
                f.blur_peak_hu_mean  = wall_stats.blur_peak_hu_mean;
                f.blur_peak_hu_std   = wall_stats.blur_peak_hu_std;
                f.blur_lung_hu_mean  = wall_stats.blur_lung_hu_mean;
                f.blur_contrast_mean = wall_stats.blur_contrast_mean;
                f.blur_contrast_std  = wall_stats.blur_contrast_std;
                f.blur_trans_width_mean = wall_stats.blur_trans_width_mean;
                f.blur_trans_width_std  = wall_stats.blur_trans_width_std;
                f.blur_edge_sharp_mean  = wall_stats.blur_edge_sharp_mean;
                f.blur_edge_sharp_std   = wall_stats.blur_edge_sharp_std;
                f.blur_contrast_gen5plus    = wall_stats.blur_contrast_gen5plus;
                f.blur_trans_width_gen5plus = wall_stats.blur_trans_width_gen5plus;
                f.blur_edge_sharp_gen5plus  = wall_stats.blur_edge_sharp_gen5plus;
                f.blur_peak_hu_gen5plus     = wall_stats.blur_peak_hu_gen5plus;
                % FWHM 版 T/D
                f.TD_fwhm_all          = wall_stats.TD_fwhm_all;
                f.TD_fwhm_std          = wall_stats.TD_fwhm_std;
                f.TD_fwhm_cv           = wall_stats.TD_fwhm_cv;
                f.TD_fwhm_gen5plus     = wall_stats.TD_fwhm_gen5plus;
                f.TD_fwhm_std_gen5plus = wall_stats.TD_fwhm_std_gen5plus;
                f.TD_fwhm_slope_vs_gen = wall_stats.TD_fwhm_slope_vs_gen;
            end
            fprintf('   管壁密度完成: mean=%.1f HU, std=%.1f\n', f.wall_hu_mean, f.wall_hu_std);
            fprintf('   边界模糊度完成: 对比=%.1f HU, 过渡带=%.2f mm, 锐度=%.2f HU/mm\n', ...
                f.blur_contrast_mean, f.blur_trans_width_mean, f.blur_edge_sharp_mean);
        catch ME
            fprintf('   管壁密度失败: %s\n', ME.message);
        end
    end

    % 写入单患者特征 CSV
    row_table = struct2table(f, 'AsArray', true);
    writetable(row_table, fullfile(FEATURES_DIR, [patient, '_airway_features.csv']));

    % 追加到汇总
    feature_table = [feature_table; row_table]; %#ok<AGROW>
    catch ME
        fprintf('  [跳过] 患者 %s 处理失败: %s\n', patient, ME.message);
        continue;
    end
end

% =========================================================================
% 4. 汇总输出
% =========================================================================
if height(feature_table) > 0
    summary_path = fullfile(FEATURES_DIR, 'airway_features_all.csv');
    writetable(feature_table, summary_path);
    fprintf('\n======================================================\n');
    fprintf('气道重塑特征计算完成！\n');
    fprintf('汇总: %s（%d 个患者）\n', summary_path, height(feature_table));
    disp(feature_table);
    fprintf('======================================================\n');
else
    fprintf('未生成任何特征（请检查 METRICS_DIR 下是否有 _full_metrics.csv）。\n');
end

% =========================================================================
% 局部函数：管壁密度/纹理 + PCA
% =========================================================================
function out = compute_wall_densitometry(ct_file, seg_file, skel_file, ...
        kernel_sizes, wall_search_mm, fwhm_bin_mm, min_tube_pts)
    % 重新加载 CT + 掩膜，构建 AQnet（优先用 AirQuant 已存骨架 skel_file，
    % 失败/缺失时退回 bwskel 自愈合），在 FWHM 测量中收集管壁 HU 统计
    out = [];

    meta_CT = niftiinfo(ct_file);
    source = double(niftiread(meta_CT));
    meta_seg = niftiinfo(seg_file);
    seg_crop = logical(niftiread(meta_seg));
    % 注意：此处不裁剪 source/seg —— ClinicalAirways 需要与完整 header (meta_CT) 尺寸匹配，
    % 与 batch_airway_quant.m 的 build_skel_and_net 保持一致；裁剪会触发构建失败。

    % 清洗掩膜
    seg_base = imfill(seg_crop, 'holes');
    CC = bwconncomp(seg_base, 26);
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, idx_max] = max(numPixels);
    seg_clean = false(size(seg_base));
    seg_clean(CC.PixelIdxList{idx_max}) = true;
    seg_base = seg_clean;

    AQnet = [];
    % 优先：使用 AirQuant 已生成、拓扑已验证的骨架（skel_output / PTKskel）
    if ~isempty(skel_file) && exist(skel_file, 'file')
        try
            meta_skel = niftiinfo(skel_file);
            skel_existing = logical(niftiread(meta_skel));
            if any(skel_existing(:))
                AQnet = ClinicalAirways(skel_existing, ...
                    'source', source, 'header', meta_CT, 'seg', seg_base, ...
                    'fillholes', 0, 'largestCC', 1, ...
                    'plane_sample_sz', 0.5, 'spline_sample_sz', 0.5);
            end
        catch
            AQnet = [];
        end
    end
    % 退回：bwskel 自愈合骨架
    if isempty(AQnet)
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
                AQnet = ClinicalAirways(skel_candidate, ...
                    'source', source, 'header', meta_CT, 'seg', seg_base, ...
                    'fillholes', 0, 'largestCC', 1, ...
                    'plane_sample_sz', 0.5, 'spline_sample_sz', 0.5);
                break;
            catch
                AQnet = [];
            end
        end
    end
    if isempty(AQnet)
        error('ClinicalAirways 构建失败');
    end

    % 标签映射
    label_matrix = zeros(size(AQnet.skel), 'uint16');
    for i = 1:length(AQnet.tubes)
        label_matrix(AQnet.tubes(i).skelpoints) = i;
    end

    aq_skel = AQnet.skel;
    aq_seg  = AQnet.seg;
    aq_source = AQnet.source;
    aq_voxdim = AQnet.voxdim;
    aq_mean_spacing = mean(aq_voxdim);
    aq_voxel_vol = prod(aq_voxdim);
    branches = aq_skel;
    dist_from_lumen = bwdist(aq_seg) * aq_mean_spacing;
    pad = ceil(wall_search_mm / aq_mean_spacing);

    % 逐分支收集：管壁 HU 统计 + 几何 + FWHM 边界模糊 + T/D（用于 PCA/异质性）
    branch_hu = struct('mean',{}, 'std',{}, 'skew',{}, 'kurt',{}, ...
                       'wt',{}, 'la',{}, 'wa',{}, ...
                       'peak_hu',{}, 'lung_hu',{}, 'contrast',{}, ...
                       'trans_w',{}, 'edge_slope',{}, ...
                       'd_in',{}, 'd_out',{}, 'td',{}, 'gen',{});
    for k = 1:length(AQnet.tubes)
        idx_pts = find(label_matrix == k);
        if length(idx_pts) < min_tube_pts
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
        other = (local_branches & (local_label ~= k) & (local_label > 0));
        D_k = bwdist(branch_k_mask);
        if any(other(:))
            D_other = bwdist(other);
            local_M = (D_k < D_other);
        else
            local_M = true(size(local_branches));
        end

        L = AQnet.tubes(k).stats.arclength;
        lumen_voxels = sum(local_M(:) & local_seg(:));
        LA = (lumen_voxels * aq_voxel_vol) / L;
        D_in = 2 * sqrt(LA / pi);

        wall_zone = local_M & ~local_seg & (local_dist <= wall_search_mm);
        hu_vals = double(local_source(wall_zone));
        d_vals  = local_dist(wall_zone);
        if isempty(hu_vals) || numel(hu_vals) < 5
            continue;
        end
        % 过滤非物理 HU（< -1024 或异常）；d_vals 必须同步过滤，
        % 否则后面 hu_vals(bin_idx==b) 的直方图逻辑索引长度失配。
        keep = (hu_vals > -1024) & (hu_vals < 1024);
        hu_vals = hu_vals(keep);
        d_vals  = d_vals(keep);
        if numel(hu_vals) < 5
            continue;
        end

        % 内径外径壁厚（FWHM 简化：用与 batch 一致的渲染值）
        edges = 0:fwhm_bin_mm:wall_search_mm;
        bin_centers = edges(1:end-1) + fwhm_bin_mm/2;
        [~, ~, bin_idx] = histcounts(d_vals, edges);
        mean_hu_bin = nan(length(bin_centers),1);
        for b = 1:length(bin_centers)
            if any(bin_idx == b)
                mean_hu_bin(b) = mean(hu_vals(bin_idx == b));
            end
        end
        mean_hu_bin = fillmissing(mean_hu_bin, 'linear');
        [peak_hu, peak_idx] = max(mean_hu_bin);
        lung_hu = min(mean_hu_bin(peak_idx:end));
        if isempty(lung_hu) || isnan(lung_hu), lung_hu = -850; end
        half_max = (peak_hu + lung_hu)/2;
        drop_idx = find(mean_hu_bin(peak_idx:end) < half_max, 1);
        if isempty(drop_idx)
            WT = NaN; trans_w = NaN; edge_slope = NaN;
        else
            WT = bin_centers(peak_idx + drop_idx - 1);
            % 外缘过渡带宽度：壁峰值 → 半高下落点的距离 (mm)
            trans_w = bin_centers(peak_idx + drop_idx - 1) - bin_centers(peak_idx);
            if trans_w > 0
                % 边界锐度 (HU/mm)：模糊/水肿 → 过渡带变宽 → 锐度减小
                edge_slope = (peak_hu - half_max) / trans_w;
            else
                edge_slope = NaN;
            end
        end
        D_out = D_in + 2*WT;
        WA_abs = pi*(D_out/2)^2 - pi*(D_in/2)^2;

        branch_hu(end+1) = struct( ...
            'mean', mean(hu_vals), ...
            'std',  std(hu_vals), ...
            'skew', skewness(hu_vals), ...
            'kurt', kurtosis(hu_vals), ...
            'wt',  WT, ...
            'la',  LA, ...
            'wa',  WA_abs, ...
            'peak_hu',  peak_hu, ...
            'lung_hu',  lung_hu, ...
            'contrast', peak_hu - lung_hu, ...
            'trans_w',  trans_w, ...
            'edge_slope', edge_slope, ...
            'd_in',  D_in, ...
            'd_out', D_out, ...
            'td',    WT / D_out, ...
            'gen',   AQnet.tubes(k).generation); %#ok<AGROW>
    end

    if isempty(branch_hu)
        error('无有效分支可统计管壁密度');
    end

    hu_mean_all = mean([branch_hu.mean]);
    hu_std_all  = mean([branch_hu.std]);
    hu_skew_all = mean([branch_hu.skew]);
    hu_kurt_all = mean([branch_hu.kurt]);
    gen_arr = [branch_hu.gen];
    out.hu_mean_all = hu_mean_all;
    out.hu_std_all  = hu_std_all;
    out.hu_skew_all = hu_skew_all;
    out.hu_kurt_all = hu_kurt_all;
    out.hu_mean_gen3 = mean([branch_hu(gen_arr==3).mean]);
    out.hu_mean_gen4 = mean([branch_hu(gen_arr==4).mean]);
    out.hu_mean_gen5 = mean([branch_hu(gen_arr==5).mean]);
    if isnan(out.hu_mean_gen3), out.hu_mean_gen3 = NaN; end
    if isnan(out.hu_mean_gen4), out.hu_mean_gen4 = NaN; end
    if isnan(out.hu_mean_gen5), out.hu_mean_gen5 = NaN; end

    % ---- FWHM 边界模糊度（急性炎性重塑：水肿/炎症 → HU 下降、对比度降低、过渡带变宽、锐度减小）----
    gen5p = gen_arr >= 5;
    out.blur_peak_hu_mean  = nanmean([branch_hu.peak_hu]);
    out.blur_peak_hu_std   = nanstd([branch_hu.peak_hu]);
    out.blur_lung_hu_mean  = nanmean([branch_hu.lung_hu]);
    out.blur_contrast_mean = nanmean([branch_hu.contrast]);
    out.blur_contrast_std  = nanstd([branch_hu.contrast]);
    out.blur_trans_width_mean = nanmean([branch_hu.trans_w]);    % 外缘过渡带 (mm)
    out.blur_trans_width_std  = nanstd([branch_hu.trans_w]);
    out.blur_edge_sharp_mean  = nanmean([branch_hu.edge_slope]); % HU/mm
    out.blur_edge_sharp_std   = nanstd([branch_hu.edge_slope]);
    % 远端小气道（>=5 代）——炎性重塑最敏感区
    out.blur_contrast_gen5plus     = nanmean([branch_hu(gen5p).contrast]);
    out.blur_trans_width_gen5plus  = nanmean([branch_hu(gen5p).trans_w]);
    out.blur_edge_sharp_gen5plus   = nanmean([branch_hu(gen5p).edge_slope]);
    out.blur_peak_hu_gen5plus      = nanmean([branch_hu(gen5p).peak_hu]);

    % ---- FWHM 版 T/D（与 CSV 分代 T/D 互相印证，反映“急剧变化”）----
    out.TD_fwhm_all  = nanmean([branch_hu.td]);
    out.TD_fwhm_std  = nanstd([branch_hu.td]);
    if out.TD_fwhm_all > 0
        out.TD_fwhm_cv = out.TD_fwhm_std / out.TD_fwhm_all;
    else
        out.TD_fwhm_cv = NaN;
    end
    out.TD_fwhm_gen5plus       = nanmean([branch_hu(gen5p).td]);
    out.TD_fwhm_std_gen5plus   = nanstd([branch_hu(gen5p).td]);
    out.TD_fwhm_slope_vs_gen   = td_slope_vs_gen(gen_arr, [branch_hu.td]);

    % PCA 异质性：每分支 [hu_mean, hu_std, hu_skew, hu_kurt, wt, la, wa, td, contrast, edge_slope]
    X = [[branch_hu.mean]', [branch_hu.std]', [branch_hu.skew]', ...
         [branch_hu.kurt]', [branch_hu.wt]', [branch_hu.la]', [branch_hu.wa]', ...
         [branch_hu.td]', [branch_hu.contrast]', [branch_hu.edge_slope]'];
    X = X(~any(isnan(X),2), :);
    if size(X,1) >= 3 && size(X,2) >= 2
        Xc = X - mean(X,1);
        Xs = std(Xc,0,1); Xs(Xs==0) = 1;
        Xn = Xc ./ Xs;
        [~, score, latent] = pca(Xn);
        out.pca_explained = latent ./ sum(latent);
        out.pca_pc1_std = std(score(:,1));
    else
        out.pca_explained = [NaN NaN NaN];
        out.pca_pc1_std = NaN;
    end
end

% =========================================================================
% 局部函数：T/D 急剧变化辅助
% =========================================================================
function s = td_slope_vs_gen(gen, tdr)
    % 每分支 T/D 对代数的稳健线性斜率（每代 T/D 变化量，反映近端→远端“急剧变化”）
    s = NaN;
    ok = ~isnan(gen) & ~isnan(tdr);
    if nnz(ok) < 3
        return;
    end
    x = double(gen(ok)); y = double(tdr(ok));
    try
        b = robustfit(x, y, 'bisquare');
        s = b(2);
    catch
        p = polyfit(x, y, 1);
        s = p(1);
    end
end

function r = td_outlier_ratio(tdr)
    % 局灶陡变：T/D 偏离中位数 |z|>2 的分支占比（MAD 稳健 z 分数）
    r = NaN;
    ok = ~isnan(tdr);
    x = tdr(ok);
    if numel(x) < 5
        return;
    end
    med = median(x);
    madv = mad(x, 1);
    if madv <= 0
        r = 0;
        return;
    end
    z = (x - med) / (1.4826 * madv);
    r = mean(abs(z) > 2);
end
