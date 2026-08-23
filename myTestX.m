%% 原生 MATLAB 气道高维组学提取脚本 (无第三方依赖版)
% 功能: 计算气管内径、外径、管壁厚度、壁面积及临床金标准 Pi10
clear; clc;

% =========================================================================
% 1. 定义专属路径 (请确保这三个文件都已存在)
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_sourceX\patient_00.nii.gz';  
seg_name   = 'D:\copd-radiomics\Airway_out\patient_00.nii_airway.nii.gz'; 
skel_name  = 'D:\copd-radiomics\PTKskel\patient_00.nii_airway.nii.gz'; 
output_csv = 'D:\copd-radiomics\Airway_out\patient_00_full_metrics.csv'; 

% =========================================================================
% 2. 加载数据与物理参数
% =========================================================================
disp('-> 正在加载 CT、掩膜与骨架数据...');
meta_CT = niftiinfo(CT_name);
source  = double(niftiread(meta_CT));
seg     = logical(niftiread(seg_name));
skel    = logical(niftiread(skel_name));

% 提取体素的物理尺寸 (毫米)
spacing = meta_CT.PixelDimensions;
mean_spacing = mean(spacing); 
voxel_vol = prod(spacing);

% =========================================================================
% 3. 气道树三维空间解析 (剥离独立分支)
% =========================================================================
disp('-> 正在进行 3D 气道分支解析...');
% 提取分支点并打断骨架，分离出独立的管道分支
bp = bwmorph3(skel, 'branchpoints');
branches = skel & ~bp;

% 标记每一个独立的分支
CC = bwconncomp(branches, 26);
label_matrix = labelmatrix(CC);

% 【避坑 2^24 报错】：仅计算距离值（第一输出），不请求索引数组
dist_from_lumen = bwdist(seg) * mean_spacing;

% =========================================================================
% 4. 核心几何组学运算 (局部 Bounding Box 法提取内外径与壁厚)
% =========================================================================
disp('-> 正在执行原生一维降维 FWHM 管壁测量 (启用自适应局部窗口)...');
results = [];
num_branches = CC.NumObjects;
pad = ceil(6 / mean_spacing); % 向外扩展 6mm 的像素个数，确保包住管壁

for k = 1:num_branches
    % 提取当前分支的骨架体素
    idx = CC.PixelIdxList{k};
    if length(idx) < 5
        continue; % 忽略过短的微小噪点分支
    end
    
    % 1. 计算当前分支的近似物理长度 (L)
    L = length(idx) * mean_spacing;
    
    % 2. 【核心修复】：局部 Bounding Box 裁切 (彻底解决 2^24 超限报错，速度提升百倍)
    [d1, d2, d3] = ind2sub(size(seg), idx);
    
    min1 = max(1, min(d1) - pad); max1 = min(size(seg,1), max(d1) + pad);
    min2 = max(1, min(d2) - pad); max2 = min(size(seg,2), max(d2) + pad);
    min3 = max(1, min(d3) - pad); max3 = min(size(seg,3), max(d3) + pad);
    
    % 裁切出微小的局部 3D 矩阵
    local_branches = branches(min1:max1, min2:max2, min3:max3);
    local_label    = label_matrix(min1:max1, min2:max2, min3:max3);
    local_seg      = seg(min1:max1, min2:max2, min3:max3);
    local_dist     = dist_from_lumen(min1:max1, min2:max2, min3:max3);
    local_source   = source(min1:max1, min2:max2, min3:max3);
    
    % 3. 划定当前分支的局部三维专属领地 (局部 Voronoi 划分)
    branch_k_mask = (local_label == k);
    other_branches_mask = (local_branches & (local_label ~= k));
    
    D_k = bwdist(branch_k_mask);
    if any(other_branches_mask(:))
        D_other = bwdist(other_branches_mask);
        local_M = (D_k < D_other); % 离自己比离别人近的区域，归属自己
    else
        local_M = true(size(local_branches)); 
    end
    
    % 4. 【内径计算】：提取属于该分支的内腔掩膜，根据体积反推内径
    lumen_voxels = sum(local_M(:) & local_seg(:));
    LumenVolume = lumen_voxels * voxel_vol;
    LA = LumenVolume / L; % 气管平均内腔面积 (Lumen Area)
    D_in = 2 * sqrt(LA / pi); % 气管平均内径 (Inner Diameter)
    Pi_val = pi * D_in; % 内部周长 (Internal Perimeter)
    
    % 5. 【管壁厚度计算】：距离-灰度统计学降维 FWHM 法
    % 提取该分支周围 5mm 内的软组织区域 (可能包含管壁)
    wall_zone = local_M & ~local_seg & (local_dist <= 5);
    d_vals = local_dist(wall_zone);
    hu_vals = local_source(wall_zone);
    
    if isempty(d_vals)
        continue;
    end
    
    % 将 3D 数据划分为 0~5mm 的同心圆刻度区间 (Binning)
    edges = 0:0.2:5;
    bin_centers = edges(1:end-1) + 0.1;
    [~, ~, bin_idx] = histcounts(d_vals, edges);
    
    mean_hu = zeros(length(bin_centers), 1);
    for b = 1:length(bin_centers)
        if any(bin_idx == b)
            mean_hu(b) = mean(hu_vals(bin_idx == b)); % 计算每个距离圈的平均 CT 灰度
        else
            mean_hu(b) = NaN;
        end
    end
    mean_hu = fillmissing(mean_hu, 'linear'); % 平滑曲线
    
    % 寻找管壁最亮处 (Peak) 和外部肺实质 (Lung Background)
    [peak_hu, peak_idx] = max(mean_hu);
    lung_hu = min(mean_hu(peak_idx:end));
    if isempty(lung_hu) || isnan(lung_hu), lung_hu = -850; end
    
    % 半高全宽阈值 (FWHM)
    half_max = (peak_hu + lung_hu) / 2;
    
    % 寻找灰度跌落至半高的位置，即为气管外壁边界
    drop_idx = find(mean_hu(peak_idx:end) < half_max, 1);
    if isempty(drop_idx)
        WT = 1.0; % 如果拟合失败，给定临床默认最小厚度
    else
        WT = bin_centers(peak_idx + drop_idx - 1); % 管壁厚度 (Wall Thickness)
    end
    
    % 6. 【外径与管壁面积计算】
    D_out = D_in + 2 * WT; % 外径
    WA = pi * (D_out/2)^2 - pi * (D_in/2)^2; % 管壁面积 (Wall Area)
    sqrt_WA = sqrt(WA);
    WA_pct = (WA / (WA + LA)) * 100; % 管壁面积占比 WA%
    
    % 保存当前分支的所有计算结果
    % [ID, 长度, 内腔面积, 壁面积, WA%, 内径, 外径, 壁厚, 内周长, 根号壁面积]
    results = [results; k, L, LA, WA, WA_pct, D_in, D_out, WT, Pi_val, sqrt_WA];
end

% =========================================================================
% 5. 导出数据表格
% =========================================================================
disp('-> 正在导出特征至 CSV...');
VariableNames = {'BranchID', 'Length_mm', 'LumenArea_mm2', 'WallArea_mm2', ...
                 'WA_pct', 'Inner_Diameter_mm', 'Outer_Diameter_mm', ...
                 'Wall_Thickness_mm', 'Pi_Perimeter_mm', 'Sqrt_WallArea'};
T = array2table(results, 'VariableNames', VariableNames);
writetable(T, output_csv);
disp(head(T, 5));
fprintf('🎉 基础几何组学提取完成！数据已保存至: %s\n', output_csv);

% =========================================================================
% 6. 计算临床金标准特征: Pi10 并自动绘制回归图
% =========================================================================
disp('-> 正在计算 COPD 黄金定量特征 (Pi10)...');

% 过滤掉异常数据 (如内径为0的噪点)
valid_idx = results(:,9) > 0 & results(:,10) > 0;
X_Pi = results(valid_idx, 9);       % X轴: 气管内周长 Pi
Y_sqrtWA = results(valid_idx, 10);  % Y轴: 管壁面积的平方根 sqrt(WA)

% 使用最小二乘法进行线性回归 (y = a*x + b)
p = polyfit(X_Pi, Y_sqrtWA, 1);

% Pi10 的数学定义: 当内周长 X = 10mm 时，理论拟合的管壁厚度参数 Y
Pi10 = polyval(p, 10);

fprintf('\n======================================================\n');
fprintf('🔥 本例患者的 Pi10 指标计算结果为: %.4f\n', Pi10);
fprintf('======================================================\n');

% 绘制专业临床回归散点图
figure('Name', 'Pi10 Regression Model', 'Color', 'w');
scatter(X_Pi, Y_sqrtWA, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'w');
hold on;
% 绘制回归直线
x_fit = linspace(min(X_Pi)*0.8, max(X_Pi)*1.2, 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% 标出 Pi=10 的黄金焦点
plot([10 10], [0 Pi10], 'k--', 'LineWidth', 1.5);
plot([0 10], [Pi10 Pi10], 'k--', 'LineWidth', 1.5);
scatter(10, Pi10, 80, 'rp', 'filled', 'MarkerEdgeColor', 'k'); % 五角星标记

% 图表美化
title(sprintf('Airway Remodeling Regression (Pi10 = %.3f)', Pi10), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Internal Perimeter (Pi, mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Square Root of Wall Area (\surdWA, mm)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2);
legend('Airway Branches', 'Regression Line', 'Pi10 Marker', 'Location', 'northwest');
hold off;

disp('🎉 全部流程执行完毕！请查看弹出的 Pi10 回归图。');