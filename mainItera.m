%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%此文件用于仿真单个主用户的情况；仿真图的横纵坐标分别为：最大可用发射功率和加权
%和速率；在仿真时，令各个次级用户权重均为1；
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%% 初始化参数
func = normalFuncSet; %导入函数集
n_monte = 10; %蒙特卡洛仿真次数
n_IRS = 5; %功率迭代的次数
acc_stop = 0.001;%主程序运行时的停止精度
max_cnt_alg = 30;%算法最大迭代次数
rate_mat = zeros(max_cnt_alg+1,n_IRS);%用于存储速率

%% 初始化场景参数
[scene,dist] = func.init();
scene.max_pow = 5;

%% 进行蒙特卡洛仿真
tic
for cnt_IRS = 1:n_IRS
    scene.m_IRS = 20*cnt_IRS;
    disp(['The number of IRS is ', num2str(scene.m_IRS),'.']);
    for cnt_monte = 1:n_monte
        disp(['Iteration = ', num2str(cnt_monte)]);
        channel = func.setChannel(scene, dist);
        %初始化反射系数矩阵和预编码矩阵
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene); 
        %计算联合信道
        [g_AP_PUs,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
        %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %%计算所有次级用户的速率和
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        rate_mat(1,cnt_IRS) = rate_mat(1,cnt_IRS) + sum_rate;
        %disp(['The front sum of rates is ',num2str(sum_rate),' bps.']);
        
        for cnt_iter = 1:max_cnt_alg
            sum_rate_tmp = sum_rate;
            %计算解码矩阵和辅助矩阵
            [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat);
            %给定其他参数的情况下计算预编码矩阵
            precode_mat = getPrecodeMat(scene,g_AP_PUs,g_AP_SUs,decode_mat,weight_mat);
            %给定其他参数的情况下计算反射系数矩阵
            reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat);
            %计算联合信道
            [g_AP_PUs,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
            %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
            [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
            %%计算所有次级用户的速率和
            sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
            rate_mat(cnt_iter+1,cnt_IRS) = rate_mat(cnt_iter+1,cnt_IRS) + sum_rate;
        end
    end
end
rate_mat = rate_mat/n_monte;
toc