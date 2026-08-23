%% 整合版：AirQuant 气道定量 (myTest4 + myTestXX + AQ1 融合)
% =========================================================================
% 功能：
%   1) 读取 CT + 气道掩膜
%   2) 自愈合骨架提取 (kernel 0/3/5/7) 构建 AQnet 拓扑树
%   3) 【修复核心】正确映射每根气管分支 (直接线性索引，替代易失败的坐标猜测雷达)
%   4) 原生 FWHM 测量：内径/外径/壁厚/管腔面积/管壁面积
%   5) 重构 3D 外壁掩膜 (label 1=内腔, label 2=外壁)
%   6) 导出 CSV (拓扑 + 几何) 与多标签 NIfTI
%   7) 计算 Pi10 并绘图
%   8) 可选肺叶分类 + 可视化 (myTest4 逻辑)
% =========================================================================
clear; clc; close all;

% =========================================================================
% 1. 路径配置
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_sourceX\patient_01.nii.gz';
seg_name   = 'D:\copd-radiomics\Airway_out\patient_01_airway.nii.gz';
skel_name  = 'D:\copd-radiomics\PTKskel\patient_01_airway.nii.gz';   % 与 seg_name 分开，避免覆盖掩膜
output_csv = 'D:\copd-radiomics\Airway_out\patient_01full_metrics.csv';
output_wall= 'D:\copd-radiomics\Airway_out\patient_01_airway_OuterWall.nii.gz';

% 测量参数
KERNEL_SIZES    = [0, 3, 5, 7];   % 自愈合断桥内核
WALL_SEARCH_MM  = 5;              % 外壁搜索最大距离 (mm)
FWHM_BIN_MM     = 0.2;            % FWHM 距离直方图 bin 宽度
MIN_TUBE_PTS    = 3;              % 少于该骨架点数则忽略该分支

% =========================================================================
% 2. 加载数据与物理参数
% =========================================================================
disp('-> 正在加载 CT 与掩膜数据...');
meta_CT = niftiinfo(CT_name);
source  = double(niftiread(meta_CT));
meta_seg = niftiinfo(seg_name);
seg_raw  = logical(niftiread(meta_seg));

spacing = meta_CT.PixelDimensions;
mean_spacing = mean(spacing);
voxel_vol = prod(spacing);

% 基础清洗：填补封闭气泡，保留最大主树
seg_base = imfill(seg_raw, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx}) = true;
seg_base = seg_clean;
%%
% =========================================================================
% 3. 自愈合骨架循环 (构建 AirQuant 拓扑树)
% =========================================================================
success = false;
for k = KERNEL_SIZES
    fprintf('\n======================================================\n');
    fprintf('🔄 正在尝试骨架净化等级: %d\n', k);

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

    disp('   -> 正在生成 3D 中心线骨架...');
    skel = bwskel(seg_for_skel, 'MinBranchLength', 15);

    disp('   -> 正在将骨架送入 AirQuant 进行图网络编译...');
    try
        AQnet = ClinicalAirways(skel, ...
            'source', source, ...
            'header', meta_CT, ...
            'seg', seg_base, ...
            'fillholes', 0, ...
            'largestCC', 1, ...
            'plane_sample_sz', 0.5, ...
            'spline_sample_sz', 0.5);

        disp('   -> 图网络编译成功！');
        niftiwrite(uint8(skel), skel_name, meta_seg, 'Compressed', true);
        success = true;
        break;
    catch ME
        fprintf('   ❌ 失败。错误信息: %s\n', ME.message);
    end
end

if ~success
    error('⚠️ 拓扑自愈合失败，掩膜存在严重粘连，请手动处理。');
end
%%
% =========================================================================
% 4. 【核心修复】分支标签映射 (在 AQnet 重定向+裁剪空间内)
% =========================================================================
% AirQuant 的 ClinicalAirways/TubeNetwork 构造时会自动对 skel/seg/source
% 做 ReorientVolume (permute+flip) 和 CropVol (裁剪到分割范围)。
% 因此 skelpoints 是「重定向+裁剪后」体积的一维线性索引
% (见 ClassifySegmentationTubes: classedskel(obj.tubes(ii).skelpoints) = ID)。
% 必须把 label_matrix 建在 AQnet.skel 上，并与 AQnet.seg / AQnet.source / AQnet.voxdim 对齐。
disp('-> 正在为每根气管分支分配标签 (在 AQnet 空间内直接线性索引映射)...');
label_matrix = zeros(size(AQnet.skel), 'uint16');
for i = 1:length(AQnet.tubes)
    label_matrix(AQnet.tubes(i).skelpoints) = i;
end
n_labeled = nnz(label_matrix > 0);
fprintf('   ✅ 已标记 %d 个骨架点，共 %d 根分支。\n', n_labeled, length(AQnet.tubes));

% =========================================================================
% 5. 原生 FWHM 几何测量 (在 AQnet 重定向+裁剪空间内)
% =========================================================================
disp('-> 正在执行 FWHM 几何测量...');

aq_skel   = AQnet.skel;      % 重定向+裁剪后的骨架
aq_seg    = AQnet.seg;       % 重定向+裁剪后的掩膜
aq_source = AQnet.source;    % 重定向+裁剪后的 CT
aq_voxdim = AQnet.voxdim;    % 重定向后的体素尺寸
aq_spacing = aq_voxdim;
aq_mean_spacing = mean(aq_voxdim);
aq_voxel_vol = prod(aq_voxdim);

branches = aq_skel;
dist_from_lumen = bwdist(aq_seg) * aq_mean_spacing;  % 到管腔边界的距离 (mm)
results = zeros(0, 9);
pad = ceil(WALL_SEARCH_MM / aq_mean_spacing);

% 外壁重构掩膜 (在 AQnet 空间)
outer_wall_mask = false(size(aq_seg));
% 各分支外壁累计，用于可视化
wall_volume_by_branch = zeros(length(AQnet.tubes), 1);

for k = 1:length(AQnet.tubes)
    idx = find(label_matrix == k);
    if length(idx) < MIN_TUBE_PTS
        continue;
    end

    [d1, d2, d3] = ind2sub(size(aq_seg), idx);
    min1 = max(1, min(d1) - pad); max1 = min(size(aq_seg,1), max(d1) + pad);
    min2 = max(1, min(d2) - pad); max2 = min(size(aq_seg,2), max(d2) + pad);
    min3 = max(1, min(d3) - pad); max3 = min(size(aq_seg,3), max(d3) + pad);

    local_branches = branches(min1:max1, min2:max2, min3:max3);
    local_label    = label_matrix(min1:max1, min2:max2, min3:max3);
    local_seg      = aq_seg(min1:max1, min2:max2, min3:max3);
    local_dist     = dist_from_lumen(min1:max1, min2:max2, min3:max3);
    local_source   = aq_source(min1:max1, min2:max2, min3:max3);

    % 局部 Voronoi 划分防粘连
    branch_k_mask = (local_label == k);
    other_branches_mask = (local_branches & (local_label ~= k) & (local_label > 0));

    D_k = bwdist(branch_k_mask);
    if any(other_branches_mask(:))
        D_other = bwdist(other_branches_mask);
        local_M = (D_k < D_other);
    else
        local_M = true(size(local_branches));
    end

    % ---- 内径 (由管腔截面积反推) ----
    L = AQnet.tubes(k).stats.arclength;
    lumen_voxels = sum(local_M(:) & local_seg(:));
    LA = (lumen_voxels * aq_voxel_vol) / L;
    D_in = 2 * sqrt(LA / pi);
    Pi_val = pi * D_in;

    % ---- 管壁 FWHM ----
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

        % 【外壁重构】：从管腔边界向外渲染到壁厚 WT 的厚度
        render_WT = max(WT, aq_mean_spacing);  % 至少渲染一层体素
        local_wall_reconstruction = local_M & ~local_seg & (local_dist <= render_WT);
        outer_wall_mask(min1:max1, min2:max2, min3:max3) = ...
            outer_wall_mask(min1:max1, min2:max2, min3:max3) | local_wall_reconstruction;
        wall_volume_by_branch(k) = nnz(local_wall_reconstruction) * aq_voxel_vol;
    end

    D_out = D_in + 2 * WT;
    WA = pi * (D_out/2)^2 - pi * (D_in/2)^2;
    sqrt_WA = sqrt(WA);
    WA_pct = (WA / (WA + LA)) * 100;

    results = [results; k, LA, WA, WA_pct, D_in, D_out, WT, Pi_val, sqrt_WA];
end

fprintf('   ✅ 完成 %d 根分支的几何测量。\n', size(results, 1));

% =========================================================================
% 6. 双表融合导出 (拓扑 + 几何)
% =========================================================================
disp('-> 正在导出组学特征 CSV...');
AQnet.ExportCSV(output_csv);
topo_table = readtable(output_csv);

VariableNames = {'ID', 'LumenArea_mm2', 'WallArea_mm2', 'WA_pct', ...
                 'Inner_Diameter_mm', 'Outer_Diameter_mm', ...
                 'Wall_Thickness_mm', 'Pi_Perimeter_mm', 'Sqrt_WallArea'};

if isempty(results)
    disp('⚠️ 警告：未能命中任何气管分支，几何厚度表为空！');
    geom_table = array2table(zeros(0, 9), 'VariableNames', VariableNames);
else
    geom_table = array2table(results, 'VariableNames', VariableNames);
end

final_table = outerjoin(topo_table, geom_table, 'Keys', 'ID', 'MergeKeys', true);
writetable(final_table, output_csv);
disp(head(final_table, 5));

% =========================================================================
% 7. 导出多标签 3D 掩膜 (1=内腔, 2=外壁) —— 映射回原始 CT 空间
% =========================================================================
disp('-> 正在把外壁掩膜从 AQnet 空间映射回原始 CT 空间...');

% --- 计算 ReorientVolume 的轴置换与翻转 (与 AirQuant 内部一致) ---
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

% --- CropVol 的裁剪范围 (AirQuant 内部用 seg 求 lims) ---
lims = AQnet.lims;  % 3x2, 每行 [min max]，在重定向后的坐标系

% --- 逆变换：把 AQnet 空间的掩膜填回完整重定向体积 ---
sz_re = size(seg_base);
sz_re = sz_re(newaxes);  % 重定向后的体积尺寸
outer_full_re = zeros(sz_re);

% CropVol 在裁剪时会先把 lims 钳制到合法范围
clims = lims;
for ii = 1:3
    if clims(ii) < 1, clims(ii) = 1; end
end
j = 1;
for ii = 4:6
    if clims(ii) > sz_re(j), clims(ii) = sz_re(j); end
    j = j + 1;
end

% 整块放置 (3D 索引)，避免逐维切片导致尺寸不匹配
r1 = clims(1,1):clims(1,2); r2 = clims(2,1):clims(2,2); r3 = clims(3,1):clims(3,2);
outer_full_re(r1, r2, r3) = outer_wall_mask;

% --- 逆翻转 ---
for i = 1:3
    if flips(i)
        outer_full_re = flip(outer_full_re, i);
    end
end

% --- 逆置换回原始空间 ---
outer_wall_orig = ipermute(outer_full_re, newaxes);

% --- 生成多标签掩膜 (原始空间) ---
outer_wall_orig = logical(outer_wall_orig);   % 转成逻辑索引
final_multi_label_mask = uint8(seg_base);              % label 1 = 内腔
final_multi_label_mask(outer_wall_orig) = 2;           % label 2 = 外壁

meta_wall = meta_seg;
meta_wall.Datatype = 'uint8';
if isfield(meta_wall, 'MultiplicativeScaling'), meta_wall.MultiplicativeScaling = 1; end
if isfield(meta_wall, 'AdditiveOffset'), meta_wall.AdditiveOffset = 0; end

niftiwrite(final_multi_label_mask, output_wall, meta_wall, 'Compressed', true);
fprintf('🎉 内外壁掩膜已保存至: %s (1=内腔, 2=外壁)\n', output_wall);

% =========================================================================
% 8. 计算临床金标准 Pi10
% =========================================================================
disp('-> 正在计算并绘制 COPD 黄金定量特征 (Pi10)...');
valid_idx = ~isnan(final_table.Pi_Perimeter_mm) & ~isnan(final_table.Sqrt_WallArea) & (final_table.Pi_Perimeter_mm > 0);
X_Pi = final_table.Pi_Perimeter_mm(valid_idx);
Y_sqrtWA = final_table.Sqrt_WallArea(valid_idx);

if length(X_Pi) > 5
    p = polyfit(X_Pi, Y_sqrtWA, 1);
    Pi10 = polyval(p, 10);

    fprintf('\n======================================================\n');
    fprintf('🔥 本例患者的 Pi10 指标计算结果为: %.4f\n', Pi10);
    fprintf('======================================================\n');

    figure('Name', 'Pi10 Regression Model', 'Color', 'w');
    scatter(X_Pi, Y_sqrtWA, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'w');
    hold on;
    x_fit = linspace(min(X_Pi)*0.8, max(X_Pi)*1.2, 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    plot([10 10], [0 Pi10], 'k--', 'LineWidth', 1.5);
    plot([0 10], [Pi10 Pi10], 'k--', 'LineWidth', 1.5);
    scatter(10, Pi10, 80, 'rp', 'filled', 'MarkerEdgeColor', 'k');
    title(sprintf('Airway Remodeling Regression (Pi10 = %.3f)', Pi10), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Internal Perimeter (Pi, mm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Square Root of Wall Area (\surdWA, mm)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    legend('Airway Branches', 'Regression Line', 'Pi10 Marker', 'Location', 'northwest');
    hold off;
else
    disp('⚠️ 有效分支数量过少，无法进行稳健的 Pi10 线性回归拟合。');
end

% =========================================================================
% 9. 肺叶分类 + 可视化 (myTest4 逻辑)
% =========================================================================
disp('-> 正在尝试肺叶分类与可视化...');
is_lobe_success = false;
try
    AQnet.ClassifyLungLobes();
    is_lobe_success = true;
    disp('   🎉 肺叶分类完成。');
catch ME
    disp('   ⚠️ 肺叶分类失败 (不影响已导出指标): ' + ME.message);
end

try
    figure('Name', 'Basic: Plot3D'); AQnet.Plot3D();
    figure('Name', 'Basic: Plot');   AQnet.Plot();
    if is_lobe_success
        figure('Name', 'Advanced: Edge weighted by generation');
        AQnet.Plot(colour='lobe', weight='generation', weightfactor=10);
    end
catch ME
    disp('   ⚠️ 可视化异常 (已跳过): ' + ME.message);
end

fprintf('\n🎉 全流程完成！\n');
fprintf('📊 CSV : %s\n', output_csv);
fprintf('🩻 掩膜: %s\n', output_wall);
