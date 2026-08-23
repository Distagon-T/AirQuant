%% 独立生成 PTKskel：仅需 CT + airway 分割，即可生成骨架 (.nii.gz)
% =========================================================================
% 用途：
%   用 MATLAB 原生 bwskel + 自愈合循环，从「气道分割掩膜」生成中心线骨架，
%   保存为与 AirQuant 兼容的 <患者名>_airway_PTKskel.nii.gz。
%   （整合版 myTest_consolidated.m 内部即用同一套逻辑，这里拆成独立脚本）
%
% 输入：
%   seg_name : 气道分割掩膜 .nii / .nii.gz（例如 Airway_out\patient_00_airway.nii.gz）
%   CT_name  : 可选。仅用于提供头信息(meta_seg 通常已足够)；若为空则直接用掩膜自身头信息。
% 输出：
%   skel_name: 生成的骨架路径（默认放在 PTKskel\ 目录）
%
% 说明：
%   PTK 官方骨架需要完整安装 Pulmonary Toolkit (tomdoel/pulmonarytoolkit)。
%   本脚本使用 MATLAB 原生 bwskel，效果等价（AirQuant 只需一个无环单连通骨架）。
%   骨架自愈合：kernel 0/3/5/7 逐级 imopen 净化，直到 ClinicalAirways 能成功编译。
% =========================================================================
clear; clc; close all;

% ---------- 1. 路径配置（按需修改） ----------
seg_name   = 'D:\copd-radiomics\Airway_out\patient_00_airway.nii.gz';   % 气道分割
CT_name    = '';   % 可选：'D:\copd-radiomics\ct_sourceX\patient_00.nii.gz'
skel_name  = 'D:\copd-radiomics\PTKskel\patient_00_airway_PTKskel.nii.gz';  % 输出骨架

% 自愈合内核大小
KERNEL_SIZES = [0, 3, 5, 7];

% ---------- 2. 加载气道分割掩膜 ----------
disp('-> 正在加载气道分割掩膜...');
assert(exist(seg_name, 'file') == 2, ['掩膜文件不存在: ', seg_name]);
meta_seg = niftiinfo(seg_name);
seg_raw  = logical(niftiread(meta_seg));

% 基础清洗：填补封闭气泡，保留最大连通主树
seg_base = imfill(seg_raw, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx}) = true;
seg_base = seg_clean;
fprintf('   ✅ 掩膜已清洗，体素数: %d\n', nnz(seg_base));

% 可选：加载 CT（用于头信息，若掩膜头信息足够可跳过）
if ~isempty(CT_name)
    meta_seg = niftiinfo(CT_name);
end

% ---------- 3. 自愈合骨架提取 ----------
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

    disp('   -> 正在生成 3D 中心线骨架 (bwskel)...');
    skel = bwskel(seg_for_skel, 'MinBranchLength', 15);

    % 验证骨架能否被 AirQuant 编译（确保无环、单连通）
    disp('   -> 正在验证骨架可被 AirQuant 编译...');
    try
        % 用一个极小占位 source 验证拓扑（无 CT 时也能测骨架）
        dummy_source = double(seg_base);
        ClinicalAirways(skel, ...
            'source', dummy_source, ...
            'seg', seg_base, ...
            'fillholes', 0, ...
            'largestCC', 1, ...
            'plane_sample_sz', 0.5, ...
            'spline_sample_sz', 0.5);
        success = true;
        break;
    catch ME
        fprintf('   ❌ 失败: %s\n', ME.message);
    end
end

if ~success
    error('⚠️ 自愈合失败，掩膜粘连严重，请先在 ITK-SNAP 中手动修整。');
end

% ---------- 4. 保存骨架 ----------
% 确保输出目录存在
[skel_dir, ~, ~] = fileparts(skel_name);
if ~isempty(skel_dir) && ~isfolder(skel_dir)
    mkdir(skel_dir);
end

meta_skel = meta_seg;
meta_skel.Datatype = 'uint8';
if isfield(meta_skel, 'MultiplicativeScaling'), meta_skel.MultiplicativeScaling = 1; end
if isfield(meta_skel, 'AdditiveOffset'), meta_skel.AdditiveOffset = 0; end

niftiwrite(uint8(skel), skel_name, meta_skel, 'Compressed', true);

fprintf('\n======================================================\n');
fprintf('🎉 PTKskel 已生成！\n');
fprintf('📄 骨架文件: %s\n', skel_name);
fprintf('🧩 骨架体素数: %d，净化等级 kernel=%d\n', nnz(skel), k);
fprintf('======================================================\n');
