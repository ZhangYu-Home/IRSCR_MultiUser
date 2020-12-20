%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%此文件用于仿真单个主用户的情况；仿真图的横纵坐标分别为：最大可用发射功率和加权
%和速率；在仿真时，令各个次级用户权重均为1；
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%% 初始化参数
func = normalFuncSet; %导入函数集
n_monte = 2; %蒙特卡洛仿真次数
n_pow = 10; %功率迭代的次数     
[scene,dist] = func.init();%初始化场景参数

%% 进行蒙特卡洛仿真
tic
for cnt_monte = 1:n_monte
    disp(['Iteration = ', num2str(cnt_monte)]);
    channel = func.setChannel(scene, dist);
    for cnt_pow = 1:n_pow
        scene.max_pow = 0.2*cnt_pow;
        %disp(['    The max Power is ', num2str(scene.max_pow),'W.']);
        %初始化反射系数矩阵和预编码矩阵
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene);       
        %计算联合信道
        [g_AP_PUs,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
        %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %计算解码矩阵和辅助矩阵
        [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat);
        %precode_mat
        
        precode_mat = getPrecodeMat(scene,g_AP_PUs,g_AP_SUs,decode_mat,weight_mat,precode_mat);
        %precode_mat
        %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %% 计算所有次级用户的速率和
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        disp(['        The sum of rates is ',num2str(sum_rate),' bps.']);
        
    end
end
toc