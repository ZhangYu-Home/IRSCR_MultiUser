%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ļ����ڷ��浥�����û������������ͼ�ĺ�������ֱ�Ϊ�������÷��书�ʺͼ�Ȩ
%�����ʣ��ڷ���ʱ��������μ��û�Ȩ�ؾ�Ϊ1��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%% ��ʼ������
func = normalFuncSet; %���뺯����
n_monte = 10; %���ؿ���������
n_pow = 10; %���ʵ����Ĵ���
acc_stop = 0.001;%����������ʱ��ֹͣ����
max_cnt_alg = 100;%�㷨����������
rate_mat = zeros(n_pow,1);%���ڴ洢����

%% ��ʼ����������
[scene,dist] = func.init();
scene.m_IRS = 80;

%% �������ؿ������
tic
for cnt_monte = 1:n_monte
    disp(['Iteration = ', num2str(cnt_monte)]);
    channel = func.setChannel(scene, dist);
    for cnt_pow = 1:n_pow
        scene.max_pow = cnt_pow;
        disp(['    The max Power is ', num2str(scene.max_pow),'W.']);
        %��ʼ������ϵ�������Ԥ�������
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene);        
        %% �̶�IRS�����ڽ����Ż��������
        %���������ŵ�
        [g_AP_PUs,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
        %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %%�������дμ��û������ʺ�
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        %disp(['The front sum of rates is ',num2str(sum_rate),' bps.']);      
        
        for cnt_iter = 1:max_cnt_alg
            sum_rate_tmp = sum_rate;
            %����������͸�������
            [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat);
            %������������������¼���Ԥ�������
            precode_mat = getPrecodeMat(scene,g_AP_PUs,g_AP_SUs,decode_mat,weight_mat);
            %������������������¼��㷴��ϵ������
            reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat);
            %���������ŵ�
            [g_AP_PUs,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
            %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
            [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
            %%�������дμ��û������ʺ�
            sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
            %disp(['The middle2 sum of rates is ',num2str(sum_rate),' bps.']);
            if(abs(sum_rate - sum_rate_tmp) < acc_stop)
                break;
            end
        end
        rate_mat(cnt_pow) = rate_mat(cnt_pow) + sum_rate;
        %disp(['The back2 sum of rates is ',num2str(sum_rate),' bps.']);
    end
end
rate_mat = rate_mat/n_monte;
save('dataPower80.mat','rate_mat');
toc