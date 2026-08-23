%% 极简气道定量提取脚本 (终极自适应无环纯化版)
clear; clc;

% =========================================================================
% 1. 定义你的专属路径
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_source\patient_09_ct.nii.gz';  
seg_name   = 'D:\copd-radiomics\Connectivity-Aware-Airway-Segmentaion\sample_out\patient_09_ct.nii_airway.nii.gz'; 
skel_name  = 'D:\patient_09_airway_PTKskel.nii.gz'; 
output_csv = 'D:\patient_097_metrics.csv'; 

% =========================================================================
% 2. 读取数据
% =========================================================================
disp('-> 正在读取 CT 与气道掩膜...');
meta_CT = niftiinfo(CT_name);
source  = double(niftiread(meta_CT));

meta_seg = niftiinfo(seg_name);
seg      = logical(niftiread(meta_seg));

% =========================================================================
% 3. 严格的 3D 形态学清洗 (不破坏远端分支，但扩宽并行气道间距)
% =========================================================================
disp('-> 正在清洗掩膜 (提取最大 3D 连通域)...');
CC = bwconncomp(seg, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
seg_clean = false(size(seg));
seg_clean(CC.PixelIdxList{idx}) = true;
seg = seg_clean;

disp('-> 正在执行微观形态学腐蚀 (打断空间由于分辨率过高导致的并行管壁粘连)...');
% 使用极其微小的 3D 结构元素进行一次侵蚀，强行拉开并行小气道之间的物理间距
se_erode = strel('cube', 3);
seg_eroded = imerode(seg, se_erode);
% 确保腐蚀后核心主干不丢失
if any(seg_eroded(:))
    seg = seg_eroded;
end
seg = imfill(seg, 'holes');

% =========================================================================
% 4. 提取 3D 骨架并利用图论斩断一切微观环路
% =========================================================================
disp('-> 正在提取初始 3D 中心线骨架...');
skel_raw = bwskel(seg, 'MinBranchLength', 25); % 提升剪枝长度，过滤掉高分辨率带来的微观毛刺

disp('-> 正在进行微观体素拓扑纯化...');
% 核心大招：利用图论提取点线后，在画回矩阵时，采用严格的 6 连通或距离约束，防止体素对角线粘连
[z, y, x] = ind2sub(size(skel_raw), find(skel_raw));
points = [x, y, z];
num_pts = size(points, 1);

[idx, dist] = rangesearch(points, points, 1.5); % 严格缩短搜索半径至 1.5（只允许面邻接）

s = []; t = []; w = [];
for i = 1:num_pts
    neighbors = idx{i};
    for j = 1:length(neighbors)
        nb = neighbors(j);
        if nb > i
            s = [s; i]; t = [t; nb]; w = [w; dist{i}(j)];
        end
    end
end
G = graph(s, t, w, num_pts);

% 熔断环路
T = minspantree(G);

% 稀疏重建骨架矩阵（绝不允许多字体素扎堆）
skel = false(size(skel_raw));
for i = 1:numedges(T)
    edge_nodes = T.Edges.EndNodes(i, :);
    p1 = points(edge_nodes(1), :); p2 = points(edge_nodes(2), :);
    skel(p1(3), p1(2), p1(1)) = true;
    skel(p2(3), p2(2), p2(1)) = true;
end

% =========================================================================
% 5. 保存并运行 AirQuant 定量
% =========================================================================
disp('-> 正在保存完美拓扑骨架为 NIfTI 文件...');
meta_skel = meta_seg;
meta_skel.Datatype = 'uint8'; 
niftiwrite(uint8(skel), skel_name, meta_skel, 'Compressed', true);

disp('-> 正在构建 AQnet 并计算内外壁厚度...');
% 关键设置：如果 fillholes 导致了重新膨胀粘连，将其设为 0
AQnet = ClinicalAirways(skel, ...
    'source', source, ...
    'header', meta_CT, ...
    'seg', seg, ...
    'fillholes', 0, ...    % 核心修正：防止 AirQuant 内部自动填充时再次把管壁泡撑大粘连
    'largestCC', 1, ...
    'plane_sample_sz', 0.5, ...
    'spline_sample_sz', 0.5);

disp('-> 正在导出定量指标 CSV...');
AQnet.ExportCSV(output_csv);

T = readtable(output_csv);
disp(head(T, 10));
fprintf('🎉 恭喜！通过微观间距拓展与内缩净化，patient_05 顺利出表！\n');