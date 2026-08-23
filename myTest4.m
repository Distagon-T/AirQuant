%%
clear; clc;

% =========================================================================
% 1. 定义路径
% =========================================================================
CT_name    = 'D:\copd-radiomics\ct_sourceX\patient_00.nii.gz';  
seg_name   = 'D:\copd-radiomics\Airway_out\patient_00.nii_airway.nii.gz'; 
skel_name  = 'D:\copd-radiomics\PTKskel\patient_00.nii_airway.nii.gz'; 
output_csv = 'D:\copd-radiomics\Airway_out\patient_00_full_metrics.csv'; 


% =========================================================================
% 
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

% AQnet.ClassifyLungLobes();


% =========================================================================
% 4. 尝试肺叶分类 (解剖学强力排雷隔离区)
% =========================================================================
disp('-> 正在尝试执行官方 AQnet.ClassifyLungLobes()...');
is_lobe_success = false;

try
    % 官方函数在遇到 patient_15 这种修剪过拓扑的病例时，内部变量会降级为 double
    % 我们用 try-catch 牢牢死守住它，允许它失败，但绝不允许它让程序崩溃！
    AQnet.ClassifyLungLobes();
    is_lobe_success = true;
    disp('   🎉 肺叶分类完全成功！将激活官方高级 RGB 可视化。');
catch ME
    disp('   ⚠️ 提示：官方肺叶分类器在当前病例内部崩溃（不影响已导出的厚度数据）。');
    fprintf('   ⚠️ 捕获到官方内部 Bug: %s\n', ME.message);
    disp('   ⚠️ 将自动屏蔽此错误并安全降级，确保后续可视化正常弹出。');
end

% =========================================================================
% 5. 可视化模块 (完全融合官方的 Basic 与 Advanced 逻辑)
% =========================================================================
disp('-> 正在生成 3D 可视化图像...');
try
    %% Basic visualisation (官方原本就不带 lobe，绝对安全，直接运行)
    % There are lots of visualisation options. Here we cover some basic options.
    figure('Name', 'Basic: Plot3D');
    AQnet.Plot3D();
    
    figure('Name', 'Basic: Plot3'); 
    AQnet.Plot3();
    
    figure('Name', 'Basic: Plot'); 
    AQnet.Plot();
    % note that the edge labels correspond to the segment indicies
    
    %% Advanced visualisation (依据分类结果，动态融合官方原有 RGB 逻辑)
    % We can use the lobe region classifications to further enrich our visualisation.
    if is_lobe_success
        % ---- 模式 A：完美状态，一字不改执行官方带有 RGB 的原代码 ----
        figure('Name', 'Advanced: Edge weighted by generation'); 
        AQnet.Plot(colour='lobe', weight='generation', weightfactor=10);
        title('Edge weighted by generation');
        
        % `Plot3D` can be resourse demanding on your specs. You may want to skip it.
        figure('Name', 'Advanced: Plot3D with Lobe'); 
        AQnet.Plot3D(colour='lobe');
        
        figure('Name', 'Advanced: PlotSpline with Lobe'); 
        AQnet.PlotSpline(colour='lobe');
    else
        % ---- 模式 B：解耦状态，剥离 colour='lobe' 的安全平替版（避开 rgb 未定义报错） ----
        figure('Name', 'Advanced (Fallback): Edge weighted by generation'); 
        AQnet.Plot(colour='generation', weight='generation', weightfactor=10);
        title('Edge weighted by generation (Fallback Mode)');
        
        figure('Name', 'Advanced (Fallback): Plot3D without Lobe'); 
        AQnet.Plot3D(); 
        
        figure('Name', 'Advanced (Fallback): PlotSpline without Lobe'); 
        AQnet.PlotSpline(); 
    end
    
catch ME
    disp(['   ⚠️ 绘图模块在显卡渲染时出现异常，已跳过部分图像。错误: ', ME.message]);
end

% =========================================================================
% 上半部分保持不变 (依然使用你那套极其稳健的自愈合掩膜与骨架提取逻辑)
% =========================================================================
% ... (从 clear; clc; 一直到 AQnet = ClinicalAirways(...) 成功构建结束) ...

% =========================================================================
% 6. 核心升级：挂载管壁几何测量模块 (提取内外径与壁厚)
% =========================================================================
disp('-> 正在启动高级几何测量引擎 (计算管壁内外径与壁厚)...');
% 注意：这一步运算量极大！算法会沿着数百条骨架向四周发射射线读取 CT 灰度
try
    % 实例化 FWHM (半高全宽) 射线测量器，传入原始 CT 和构建好的网络树
    measurer = measure.AirwayFWHMesl('source', source, ...
        'header', meta_CT, ...
        'net', AQnet);

    disp('   -> 射线发射与边界拟合中 (可能需要几分钟)...');
    measurer.Measure(); 
    disp('   🎉 壁厚测量完成！');

    % 将计算出的厚度数据强行绑定回网络节点中
    AQnet.measurements = measurer;
    is_measure_success = true;
catch ME
    warning('射线测量引擎崩溃，可能由于内存不足或管腔太细。跳过厚度计算。错误信息: %s', ME.message);
    is_measure_success = false;
end




% =========================================================================
% 6. 最终结束汇报
% =========================================================================
fprintf('\n======================================================\n');
disp('-> 正在重新验证定量指标 CSV...');
if exist(output_csv, 'file')
    T = readtable(output_csv);
    disp(head(T, 5)); % 打印前 5 行确保数据完整
    fprintf('数据安全保存在: %s\n', output_csv);
else
    warning('未检测到 CSV 文件，请检查上文循环中 ExportCSV 是否报错。');
end