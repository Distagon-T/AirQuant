%% 高级气道形态与拓扑组学特征提取 (Advanced Airway Phenotyping) — 融合版
% =========================================================================
% 依赖：先运行 batch_airway_quant.m，确保已生成 _full_metrics.csv 和
%       _airway_OuterWall.nii.gz（label 1=内腔, 2=外壁）。
%
% 三大类特征：
%   一、形态学分级：各代 WA% 均值、T/D (WT/OD)、内径/外径、分支数
%   二、图论拓扑：分叉角(parent/sibling)、迂曲度、终端分支/修剪率、最大代数
%   三、管壁密度纹理：复用 OuterWall 掩膜(==2)提取 HU mean/std/skew/kurt+分位数
%       （无需重建 AQnet，快且省内存），并对队列特征做 PCA 异质性降维
%
% 输出：
%   - <FEATURES_DIR>/<患者>_airway_features.csv       （每患者一行）
%   - <FEATURES_DIR>/airway_features_all.csv          （汇总矩阵）
%   - <METRICS_DIR>/COPD_Cohort_Advanced_Phenotypes.csv （你指定的队列总表）
%   - <FEATURES_DIR>/pca_report.txt                   （PCA 解释方差报告）
% =========================================================================
clear; clc; close all;

% =========================================================================
% 1. 路径配置
% =========================================================================
CT_INPUT_DIR = 'E:\DICOM\2026-04-niftitmp';          % 原始 CT 目录 (用于提取 HU 密度)
METRICS_DIR  = 'E:\DICOM\2026-04-Airway_metrics';    % batch_airway_quant 输出目录
FEATURES_DIR = 'E:\DICOM\2026-04-Airway_features';   % 本脚本特征输出目录
OUTPUT_CSV   = fullfile(METRICS_DIR, 'COPD_Cohort_Advanced_Phenotypes.csv');

mkdir(FEATURES_DIR);

% =========================================================================
% 2. 收集所有处理成功的患者文件夹
% =========================================================================
d = dir(fullfile(METRICS_DIR, '*_airquant'));
patient_folders = {d.name};
patient_folders = patient_folders(cellfun(@(x) isdir(fullfile(METRICS_DIR,x)), patient_folders));
fprintf('发现 %d 个患者数据包，开始提取高级组学特征...\n', numel(patient_folders));

cohort_features = table();

% =========================================================================
% 3. 逐个患者提取特征
% =========================================================================
for i = 1:numel(patient_folders)
    patient_dir  = patient_folders{i};
    patient_name = patient_dir(1:end-9);   % 去掉末尾 '_airquant' (更稳妥，避免 strrep 误伤)
    fprintf('[%d/%d] 分析患者: %s ...\n', i, numel(patient_folders), patient_name);

    csv_file  = fullfile(METRICS_DIR, patient_dir, [patient_name, '_full_metrics.csv']);
    wall_file = fullfile(METRICS_DIR, patient_dir, [patient_name, '_airway_OuterWall.nii.gz']);
    info_file = fullfile(METRICS_DIR, patient_dir, [patient_name, '_airquant_info.json']);

    % --- 从 info JSON 读 CT 路径 (确保用对序列, 而非目录第一个) 和 Pi10 ---
    ct_file = ''; Pi10 = NaN;
    if exist(info_file, 'file')
        info = jsondecode(fileread(info_file));
        if isfield(info, 'selected_ct'), ct_file = info.selected_ct; end
        if isfield(info, 'Pi10'),        Pi10    = info.Pi10;        end
    end
    if isempty(ct_file) || ~exist(ct_file, 'file')
        ct_search = dir(fullfile(CT_INPUT_DIR, patient_name, '*.nii.gz'));
        if ~isempty(ct_search)
            ct_file = fullfile(ct_search(1).folder, ct_search(1).name);
        end
    end

    % 检查文件完整性
    if ~exist(csv_file, 'file') || ~exist(wall_file, 'file')
        fprintf('  ⚠️ 缺少指标表或外壁掩膜，跳过。\n');
        continue;
    end

    % ---------------------------------------------------------------------
    % A. 读取数据
    % ---------------------------------------------------------------------
    T = readtable(csv_file);
    valid_idx = ~isnan(T.Inner_Diameter_mm) & (T.Inner_Diameter_mm > 0);
    T_valid = T(valid_idx, :);
    if isempty(T_valid)
        continue;
    end

    % 确保必要列存在
    need = {'generation','WA_pct','Wall_Thickness_mm','Outer_Diameter_mm', ...
            'Inner_Diameter_mm','stats_tortuosity','stats_sibling_deg','stats_parent_deg'};
    for n = need
        if ~ismember(n{1}, T_valid.Properties.VariableNames)
            T_valid.(n{1}) = nan(height(T_valid), 1);
        end
    end

    % ---------------------------------------------------------------------
    % B. 第一类：形态学分级特征
    % ---------------------------------------------------------------------
    T_valid.T_D_ratio = T_valid.Wall_Thickness_mm ./ T_valid.Outer_Diameter_mm;

    gen = T_valid.generation;
    gen_3_idx = (gen == 3);
    gen_4_idx = (gen == 4);
    gen_5_idx = (gen >= 5);    % 5代及以上的远端小气道

    WA_pct_Gen3     = nanmean(T_valid.WA_pct(gen_3_idx));
    WA_pct_Gen4     = nanmean(T_valid.WA_pct(gen_4_idx));
    WA_pct_Gen5plus = nanmean(T_valid.WA_pct(gen_5_idx));
    WA_pct_all      = nanmean(T_valid.WA_pct);

    TD_ratio_Gen4    = nanmean(T_valid.T_D_ratio(gen_4_idx));
    TD_ratio_Gen5plus= nanmean(T_valid.T_D_ratio(gen_5_idx));
    TD_ratio_all     = nanmean(T_valid.T_D_ratio);

    Din_mean = nanmean(T_valid.Inner_Diameter_mm);
    Dout_mean= nanmean(T_valid.Outer_Diameter_mm);

    % ---------------------------------------------------------------------
    % C. 第二类：图论与拓扑网络特征
    % ---------------------------------------------------------------------
    Branch_Count_Gen4    = sum(gen_4_idx);
    Branch_Count_Gen5plus= sum(gen_5_idx);

    Mean_Tortuosity     = nanmean(T_valid.stats_tortuosity);
    Mean_Sibling_Angle  = nanmean(T_valid.stats_sibling_deg);
    Mean_Parent_Angle   = nanmean(T_valid.stats_parent_deg);

    % 终端分支（无子代）
    hasch = false(height(T_valid), 1);
    if ismember('children_1', T_valid.Properties.VariableNames)
        hasch = hasch | ~isnan(T_valid.children_1);
    end
    if ismember('children_2', T_valid.Properties.VariableNames)
        hasch = hasch | ~isnan(T_valid.children_2);
    end
    terminal = ~hasch;
    n_terminal_total = sum(terminal);
    n_terminal_gen5plus = sum((gen >= 5) & terminal);
    pruning_ratio_gen5 = n_terminal_gen5plus / max(1, height(T_valid));
    max_generation = max(gen);

    % ---------------------------------------------------------------------
    % D. 第三类：管壁微观密度与纹理异质性 (复用 OuterWall 掩膜, 快且省内存)
    % ---------------------------------------------------------------------
    Wall_Mean_HU = NaN; Wall_Std_HU = NaN; Wall_Skew = NaN; Wall_Kurt = NaN;
    Wall_P5 = NaN; Wall_P25 = NaN; Wall_P50 = NaN; Wall_P75 = NaN; Wall_P95 = NaN;

    if ~isempty(ct_file) && exist(ct_file, 'file')
        try
            ct_vol    = double(niftiread(niftiinfo(ct_file)));
            wall_mask = niftiread(niftiinfo(wall_file));   % 1=内腔, 2=外壁
            if ~isequal(size(ct_vol), size(wall_mask))
                warning('CT 与外壁掩膜尺寸不一致，跳过 HU。');
            else
                wall_voxels_hu = ct_vol(wall_mask == 2);
                wall_voxels_hu = wall_voxels_hu(wall_voxels_hu > -1024 & wall_voxels_hu < 1024);
                if ~isempty(wall_voxels_hu)
                    Wall_Mean_HU = mean(wall_voxels_hu);
                    Wall_Std_HU  = std(wall_voxels_hu);
                    Wall_Skew    = skewness(wall_voxels_hu);
                    Wall_Kurt    = kurtosis(wall_voxels_hu);
                    prct = prctile(wall_voxels_hu, [5 25 50 75 95]);
                    Wall_P5 = prct(1); Wall_P25 = prct(2); Wall_P50 = prct(3);
                    Wall_P75 = prct(4); Wall_P95 = prct(5);
                end
            end
        catch ME
            fprintf('  ⚠️ 读取 %s 的 CT/掩膜计算 HU 失败: %s\n', patient_name, ME.message);
        end
    end

    % ---------------------------------------------------------------------
    % E. 构建单病历特征行
    % ---------------------------------------------------------------------
    patient_row = table( ...
        {patient_name}, height(T_valid), Pi10, ...
        WA_pct_Gen3, WA_pct_Gen4, WA_pct_Gen5plus, WA_pct_all, ...
        TD_ratio_Gen4, TD_ratio_Gen5plus, TD_ratio_all, ...
        Din_mean, Dout_mean, ...
        Branch_Count_Gen4, Branch_Count_Gen5plus, ...
        n_terminal_total, n_terminal_gen5plus, pruning_ratio_gen5, max_generation, ...
        Mean_Tortuosity, Mean_Sibling_Angle, Mean_Parent_Angle, ...
        Wall_Mean_HU, Wall_Std_HU, Wall_Skew, Wall_Kurt, ...
        Wall_P5, Wall_P25, Wall_P50, Wall_P75, Wall_P95, ...
        'VariableNames', { ...
        'Patient_ID', 'Valid_Branches_Count', 'Pi10', ...
        'WA_pct_Gen3', 'WA_pct_Gen4', 'WA_pct_Gen5plus', 'WA_pct_all', ...
        'TD_ratio_Gen4', 'TD_ratio_Gen5plus', 'TD_ratio_all', ...
        'Din_mean', 'Dout_mean', ...
        'Branch_Count_Gen4', 'Branch_Count_Gen5plus', ...
        'N_Terminal_Total', 'N_Terminal_Gen5plus', 'Pruning_Ratio_Gen5', 'Max_Generation', ...
        'Mean_Tortuosity', 'Mean_Sibling_Angle', 'Mean_Parent_Angle', ...
        'Wall_Mean_HU', 'Wall_Std_HU', 'Wall_Skew_HU', 'Wall_Kurt_HU', ...
        'Wall_HU_P5', 'Wall_HU_P25', 'Wall_HU_P50', 'Wall_HU_P75', 'Wall_HU_P95'});

    cohort_features = [cohort_features; patient_row]; %#ok<AGROW>

    % 每患者单独特征文件
    writetable(patient_row, fullfile(FEATURES_DIR, [patient_name, '_airway_features.csv']));
end

% =========================================================================
% 4. 队列汇总输出
% =========================================================================
if isempty(cohort_features)
    disp('⚠️ 没有成功提取到任何特征。');
    return;
end

% 汇总矩阵
all_csv = fullfile(FEATURES_DIR, 'airway_features_all.csv');
writetable(cohort_features, all_csv);

% 用户指定的队列总表
writetable(cohort_features, OUTPUT_CSV);

% ---------------------------------------------------------------------
% 5. PCA 异质性降维（跨患者，对数值特征矩阵）
% ---------------------------------------------------------------------
feat_names = {'WA_pct_Gen4','WA_pct_Gen5plus','TD_ratio_Gen4','TD_ratio_Gen5plus', ...
              'Mean_Tortuosity','Mean_Sibling_Angle','Mean_Parent_Angle', ...
              'Wall_Mean_HU','Wall_Std_HU','Wall_Skew_HU','Wall_Kurt_HU'};
feat_present = feat_names(ismember(feat_names, cohort_features.Properties.VariableNames));
if numel(patient_folders) >= 3 && ~isempty(feat_present)
    X = table2array(cohort_features(:, feat_present));
    % 仅用有限值行做标准化+PCA（避免 NaN 污染）
    ok = all(~isnan(X), 2);
    if sum(ok) >= 3
        Xs = X(ok, :);
        Xc = Xs - mean(Xs, 1);
        sd = std(Xc, 0, 1); sd(sd == 0) = 1;
        Xn = Xc ./ sd;
        [coeff, score, latent] = pca(Xn);
        explained = latent ./ sum(latent);

        % 报告
        fid = fopen(fullfile(FEATURES_DIR, 'pca_report.txt'), 'w', 'n', 'UTF-8');
        fprintf(fid, 'Airway Phenotype PCA Report\n');
        fprintf(fid, '==========================\n');
        fprintf(fid, 'Features used: %s\n', strjoin(feat_present, ', '));
        fprintf(fid, 'Patients used: %d / %d\n\n', sum(ok), height(cohort_features));
        fprintf(fid, 'Explained variance ratio:\n');
        for c = 1:min(5, numel(explained))
            fprintf(fid, '  PC%d: %.4f (cum %.4f)\n', c, explained(c), sum(explained(1:c)));
        end
        fprintf(fid, '\nPC1 loadings:\n');
        for c = 1:numel(feat_present)
            fprintf(fid, '  %-24s %.4f\n', feat_present{c}, coeff(c,1));
        end
        fclose(fid);

        % 把前 3 个主成分附加到汇总表
        score_all = nan(height(cohort_features), 3);
        score_all(ok, :) = score(:, 1:3);
        cohort_features.PC1 = score_all(:,1);
        cohort_features.PC2 = score_all(:,2);
        cohort_features.PC3 = score_all(:,3);
        writetable(cohort_features, all_csv);
        writetable(cohort_features, OUTPUT_CSV);
        fprintf('\nPCA 完成。前3主成分解释方差: %.2f%% / %.2f%% / %.2f%%\n', ...
            explained(1)*100, explained(2)*100, explained(3)*100);
    else
        fprintf('\n⚠️ 有效患者数不足，跳过 PCA。\n');
    end
else
    fprintf('\n⚠️ 患者数不足或特征缺失，跳过 PCA。\n');
end

% =========================================================================
% 6. 完成
% =========================================================================
fprintf('\n======================================================\n');
fprintf('提取完成！队列高级特征总表已保存至:\n%s\n', OUTPUT_CSV);
fprintf('汇总矩阵: %s\n', all_csv);
fprintf('PCA 报告: %s\n', fullfile(FEATURES_DIR, 'pca_report.txt'));
disp(head(cohort_features, min(5, height(cohort_features))));
fprintf('======================================================\n');
