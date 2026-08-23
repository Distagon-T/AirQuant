
%% 极简气道定量提取脚本 (彻底修复 BFS 拓扑环路报错版)
clear; clc;

% =========================================================================
% 1. 定义你的专属路径
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_source\patient_05_ct.nii.gz';  
seg_name   = 'D:\copd-radiomics\Connectivity-Aware-Airway-Segmentaion\sample_out\patient_05_ct.nii_airway.nii.gz'; 
skel_name  = 'D:\patient_05_airway_PTKskel.nii.gz'; 
output_csv = 'D:\patient_05_metrics.csv'; 

% =========================================================================
% 2. 读取数据
% =========================================================================
disp('-> 正在读取 CT 与气道掩膜...');
meta_CT = niftiinfo(CT_name);
source  = double(niftiread(meta_CT));

meta_seg = niftiinfo(seg_name);
seg      = logical(niftiread(meta_seg));

% =========================================================================
% 3. 终极 3D 形态学清洗 (破除"甜甜圈"环路陷阱)
% =========================================================================
disp('-> 正在清洗掩膜 (提取最大 3D 连通域)...');
CC = bwconncomp(seg, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
seg_clean = false(size(seg));
seg_clean(CC.PixelIdxList{idx}) = true;
seg = seg_clean;

% 【核心修复 A】：3D 形态学闭运算（糊上微小隧道，强行打断环路）
disp('-> 正在进行 3D 形态学闭运算 (修补管壁微小隧道，防止骨架打结)...');
% 使用 3x3x3 的立方体结构元素，填补所有细小缝隙
se_close = strel('cube', 3); 
seg = imclose(seg, se_close);

% 【核心修复 B】：3D 形态学开运算（打磨表面毛刺）
disp('-> 正在进行 3D 形态学开运算 (去除表面毛刺和假阳性凸起)...');
se_open = strel('sphere', 1);
seg = imopen(seg, se_open);

disp('-> 正在填补全封闭 3D 管壁空洞...');
seg = imfill(seg, 'holes');

% =========================================================================
% 4. 提取并清理 3D 骨架
% =========================================================================
disp('-> 正在生成 3D 中心线骨架 (剔除过短的假树枝)...');
% 提高 MinBranchLength，强行剪掉短于 20 个体素的冗余毛刺分支
skel = bwskel(seg, 'MinBranchLength', 20);

% 确保骨架单一连通
CC_skel = bwconncomp(skel, 26);
numPixels_skel = cellfun(@numel, CC_skel.PixelIdxList);
[~, idx_skel] = max(numPixels_skel);
skel_clean = false(size(skel));
skel_clean(CC_skel.PixelIdxList{idx_skel}) = true;
skel = skel_clean;

% =========================================================================
% 5. 保存与执行 AirQuant
% =========================================================================
disp('-> 正在保存骨架为 NIfTI 文件...');
meta_skel = meta_seg;
meta_skel.Datatype = 'uint8'; 
niftiwrite(uint8(skel), skel_name, meta_skel, 'Compressed', true);

disp('-> 正在构建 AQnet 并计算内外壁厚度 (请耐心等待)...');
% 如果这里再报错，说明你需要使用 AirQuant 自带的 Skeleton3D 函数
AQnet = ClinicalAirways(skel, ...
    'source', source, ...
    'header', meta_CT, ...
    'seg', seg, ...
    'fillholes', 1, ...
    'largestCC', 1, ...
    'plane_sample_sz', 0.5, ...
    'spline_sample_sz', 0.5);

AQnet.ClassifyLungLobes();

%% Basic visualisation
% There are lots of visualisation options. Here we cover some basic options.

figure;
AQnet.Plot3D();

figure; AQnet.Plot3();
figure; AQnet.Plot();
% note that the edge labels correspond to the segment indicies

%% Advanced visualisation
% We can use the lobe region classifications to further enrich our visualisation.
% Checkout the documentation for more advanced visualision use cases.
figure; AQnet.Plot(colour='lobe', weight='generation', weightfactor=10);
title('Edge weighted by generation')
figure; AQnet.Plot3D(colour='lobe');
% `Plot3D` can be resourse demanding on your specs. You maywant to
% skip it.
figure; AQnet.PlotSpline(colour='lobe');
disp('-> 正在导出定量指标 CSV...');
AQnet.ExportCSV(output_csv);

T = readtable(output_csv);
disp(head(T, 10));
disp(['🎉 恭喜！全部气道定量数据提取完成，已保存至: ', output_csv]);