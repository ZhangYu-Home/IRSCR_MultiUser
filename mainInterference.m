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
acc_stop = 0.01;%����������ʱ��ֹͣ����
max_cnt_alg = 100;%�㷨����������
rate_mat = zeros(n_pow,3);%���ڴ洢����

%% ��ʼ����������
[scene,dist] = func.init();
scene.m_IRS = 40;
scene.max_pow = 5;

%% �������ؿ������
tic
%����ʱ���������������
rand('seed',sum(100*clock));
for cnt_monte = 1:n_monte
    disp(['Iteration = ', num2str(cnt_monte)]);
    channel = func.setChannel(scene, dist);
    for cnt_pow = 1:n_pow
        %scene.leak_pow = 2*10^(cnt_pow-5);
        scene.leak_pow = (1e-5) * cnt_pow;
        disp(['The leak Power is ', num2str(scene.leak_pow),'W.']);
        %��ʼ������ϵ�������Ԥ�������
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene); 
         
        %% ��IRS�����ڽ����Ż��������
        %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
        [sig_mat,jam_mat] = func.getSigAndJamMat(channel.h_AP_SUs,precode_mat,scene.noise_SU);
        %%�������дμ��û������ʺ�
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        %disp(['The front sum of rates is ',num2str(sum_rate),' bps.']);      
        for cnt_iter = 1:max_cnt_alg
            sum_rate_tmp = sum_rate;
            %����������͸�������
            [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,channel.h_AP_SUs,precode_mat);
            %������������������¼���Ԥ�������
            precode_mat = getPrecodeMat(scene,channel.h_AP_PUs,channel.h_AP_SUs,decode_mat,weight_mat);
            %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
            [sig_mat,jam_mat] = func.getSigAndJamMat(channel.h_AP_SUs,precode_mat,scene.noise_SU);
            %�������дμ��û������ʺ�
            sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
            %disp(['The middle1 sum of rates is ',num2str(sum_rate),' bps.']);
            if(abs(sum_rate - sum_rate_tmp) < acc_stop)
                break;
            end
        end
        rate_mat(cnt_pow,1) = rate_mat(cnt_pow,1) + sum_rate;
        %disp(['The back1 sum of rates is ',num2str(sum_rate),' bps.']);
        
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
            %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
            [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
            %�������дμ��û������ʺ�
            sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
            %disp(['The middle2 sum of rates is ',num2str(sum_rate),' bps.']);
            if(abs(sum_rate - sum_rate_tmp) < acc_stop)
                break;
            end
        end
        rate_mat(cnt_pow,2) = rate_mat(cnt_pow,2) + sum_rate;
        %disp(['The back2 sum of rates is ',num2str(sum_rate),' bps.']);
        
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
            %disp(['The middle3 sum of rates is ',num2str(sum_rate),' bps.']);
            if(abs(sum_rate - sum_rate_tmp) < acc_stop)
                break;
            end
        end
        rate_mat(cnt_pow,3) = rate_mat(cnt_pow,3) + sum_rate;
        %disp(['The back3 sum of rates is ',num2str(sum_rate),' bps.']);
    end
end
rate_mat = rate_mat/n_monte;
save('dataInterference.mat','rate_mat');
toc