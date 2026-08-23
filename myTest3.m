%% 极简气道定量提取脚本 (工业级 Try-Catch 自愈合防弹版)
clear; clc;

% =========================================================================
% 1. 定义你的专属路径
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_sourceX\patient_01_ct.nii.gz';  
seg_name   = 'D:\copd-radiomics\Airway_out\patient_01_ct.nii_airway.nii.gz'; 
skel_name  = 'D:\copd-radiomics\PTKskel\patient_01_airway_PTKskel.nii.gz'; 
output_csv = 'D:\copd-radiomics\Airway_out\patient_01_metrics.csv'; 
% =========================================================================
% 2. 读取数据 (打好基础)
% =========================================================================
disp('-> 正在读取 CT 与气道掩膜...');
meta_CT = niftiinfo(CT_name);
source  = double(niftiread(meta_CT));

meta_seg = niftiinfo(seg_name);
seg_raw  = logical(niftiread(meta_seg));

% 基础清洗：填补绝对封闭的内部气泡，保留最大主树
seg_base = imfill(seg_raw, 'holes');
CC = bwconncomp(seg_base, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
seg_clean = false(size(seg_base));
seg_clean(CC.PixelIdxList{idx}) = true;
seg_base = seg_clean;

% =========================================================================
% 3. 开启自愈合循环 (Auto-Healing Loop)
% =========================================================================
% 设置逐渐增强的"断桥"内核大小 (0 表示原图，3/5/7 表示逐渐加强的切断力)
kernel_sizes = [0, 3, 5, 7]; 
success = false;

for k = kernel_sizes
    fprintf('\n======================================================\n');
    fprintf('🔄 正在尝试骨架净化等级: %d (内核大小)\n', k);
    
    % --- 步骤 A：为提取骨架制作"替身掩膜" ---
    if k == 0
        seg_for_skel = seg_base; % 第一次尝试不破坏任何结构
    else
        % 使用开运算(Open)：它能精准斩断细小的假桥，但完全保留粗大的主气管(不伤隆突)
        se = strel('cube', k);
        seg_for_skel = imopen(seg_base, se);
        
        % 重新锁定最大的主干，丢弃被斩断的孤岛
        CC_skel = bwconncomp(seg_for_skel, 26);
        numPixels_skel = cellfun(@numel, CC_skel.PixelIdxList);
        [~, idx_skel] = max(numPixels_skel);
        temp = false(size(seg_for_skel));
        temp(CC_skel.PixelIdxList{idx_skel}) = true;
        seg_for_skel = temp;
    end
    
    % --- 步骤 B：生成 3D 骨架 ---
    disp('   -> 正在生成 3D 中心线骨架...');
    skel = bwskel(seg_for_skel, 'MinBranchLength', 15);
    
    % --- 步骤 C：尝试运行 AirQuant (核心排雷区) ---
    disp('   -> 正在将骨架送入 AirQuant 进行图网络编译...');
    try
        % 【关键】：骨架用的是净化后的 skel，但边界和管壁厚度依然用最原始的 seg_base！
        AQnet = ClinicalAirways(skel, ...
            'source', source, ...
            'header', meta_CT, ...
            'seg', seg_base, ... 
            'fillholes', 0, ...  % 禁止 AirQuant 乱动我们的原始掩膜
            'largestCC', 1, ...
            'plane_sample_sz', 0.5, ...
            'spline_sample_sz', 0.5);
            
        disp('   -> 图网络编译成功！正在导出定量指标 CSV...');
        AQnet.ExportCSV(output_csv);
        success = true;
        
        % 保存成功的完美骨架供以后查看
        meta_skel = meta_seg;
        meta_skel.Datatype = 'uint8'; 
        niftiwrite(uint8(skel), skel_name, meta_skel, 'Compressed', true);
        
        break; % 成功出表，跳出循环！
        
    catch ME
        % 如果报错 (比如经典的 BFS not expected)，捕获错误但不崩溃！
        fprintf('   ❌ 失败：当前等级无法消除图论环路。错误信息: %s\n', ME.message);
        disp('   -> 自动启动防弹机制，准备进入下一级更强的净化处理...');
    end
end

% =========================================================================
% 4. 结果汇报
% =========================================================================
fprintf('\n======================================================\n');
if success
    T = readtable(output_csv);
    disp(head(T, 10));
    fprintf('🎉 破局成功！自愈合管线完美跑通，指标已安全导出至: %s\n', output_csv);
else
    fprintf('⚠️ 经过所有等级的尝试，依然存在极端拓扑错误。请在 ITK-SNAP 中手动切断明显的粘连。\n');
end

AQnet.ClassifyLungLobes();


% =========================================================================
% 4. 尝试肺叶分类 (解剖学 RGB 染色的前置条件)
% =========================================================================
disp('-> 正在尝试执行官方 AQnet.ClassifyLungLobes()...');
is_lobe_success = false;
try
    AQnet.ClassifyLungLobes();
    is_lobe_success = true;
    disp('   🎉 肺叶分类完全成功！将激活官方高级 RGB 可视化。');
catch
    disp('   ⚠️ 提示：由于拓扑被自愈合管线修剪，标准肺叶分类未能完全匹配。');
    disp('   ⚠️ 将自动降级为安全染色模式，不影响厚度数据导出。');
end

% =========================================================================
% 5. 可视化模块 (完全融合官方的 Basic 与 Advanced 逻辑)
% =========================================================================
disp('-> 正在生成 3D 可视化图像...');
try
    %% Basic visualisation (这部分官方原本就不带 lobe，绝对安全，直接运行)
    % There are lots of visualisation options. Here we cover some basic options.
    figure('Name', 'Basic: Plot3D');
    AQnet.Plot3D();

    figure('Name', 'Basic: Plot3'); 
    AQnet.Plot3();

    figure('Name', 'Basic: Plot'); 
    AQnet.Plot();
    % note that the edge labels correspond to the segment indicies

    %% Advanced visualisation (根据是否成功分类，动态执行官方代码)
    % We can use the lobe region classifications to further enrich our visualisation.
    if is_lobe_success
        % ---- 完美状态：一字不改执行官方原代码 ----
        figure('Name', 'Advanced: Edge weighted by generation'); 
        AQnet.Plot(colour='lobe', weight='generation', weightfactor=10);
        title('Edge weighted by generation');

        % `Plot3D` can be resourse demanding on your specs. You may want to skip it.
        figure('Name', 'Advanced: Plot3D with Lobe'); 
        AQnet.Plot3D(colour='lobe');

        figure('Name', 'Advanced: PlotSpline with Lobe'); 
        AQnet.PlotSpline(colour='lobe');

    else
        % ---- 降级状态：剥离 colour='lobe' 的安全平替版 ----
        figure('Name', 'Advanced (Fallback): Edge weighted by generation'); 
        AQnet.Plot(colour='generation', weight='generation', weightfactor=10);
        title('Edge weighted by generation (Fallback Mode)');

        % 降级为不带颜色的基础 3D 和样条图
        figure('Name', 'Advanced (Fallback): Plot3D'); 
        AQnet.Plot3D(); 

        figure('Name', 'Advanced (Fallback): PlotSpline'); 
        AQnet.PlotSpline(); 
    end

catch ME
    disp(['   ⚠️ 绘图模块在渲染时出现异常，已跳过。错误信息: ', ME.message]);
end

% =========================================================================
% 6. 最终数据导出
% =========================================================================
disp('-> 正在导出定量指标 CSV...');
AQnet.ExportCSV(output_csv);

T = readtable(output_csv);
disp(head(T, 10));
disp(['🎉 恭喜！数据提取与自适应可视化完成，已保存至: ', output_csv]);