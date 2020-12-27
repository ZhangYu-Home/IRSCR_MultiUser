%% 计算IRS反射系数
function reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat)
    % 计算所有预编码矩阵与其共轭转职乘积的和
    Q_s = zeros(scene.n_ante_AP,scene.n_ante_AP);
    for i = 1:scene.n_SU
        Q_s = Q_s + precode_mat(:,:,i)*precode_mat(:,:,i)';
    end
    
    % 计算次级用户相关的参数
    C = channel.h_AP_IRS*Q_s*channel.h_AP_IRS';
    B_0 = zeros(scene.m_IRS,scene.m_IRS);
    D_0 = zeros(scene.m_IRS,scene.m_IRS);
    for i = 1:scene.n_SU
        B_0 = B_0 + channel.h_IRS_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 + channel.h_AP_IRS*Q_s*channel.h_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 - channel.h_AP_IRS*precode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
    end
    d_0_conj = conj(diag(D_0));
    gamma_0 = B_0.*(C.');
    
    %计算主用户相关的参数
    B_k = zeros(scene.m_IRS,scene.m_IRS,scene.n_PU);
    D_k = zeros(scene.m_IRS,scene.m_IRS,scene.n_PU);
    leak_pow_k = zeros(scene.n_PU,1);
    d_k_conj = zeros(scene.m_IRS,scene.n_PU);
    gamma_k = zeros(scene.m_IRS,scene.m_IRS,scene.n_PU);
    for i = 1:scene.n_PU
        B_k(:,:,i) = channel.h_IRS_PUs(:,:,i)'*channel.h_IRS_PUs(:,:,i);
        D_k(:,:,i) = channel.h_AP_IRS*Q_s*channel.h_AP_PUs(:,:,i)'*channel.h_IRS_PUs(:,:,i);
        leak_pow_k(i) = scene.leak_pow-real(trace(channel.h_AP_PUs(:,:,i)*Q_s*channel.h_AP_PUs(:,:,i)'));
        d_k_conj(:,i) = conj(diag(D_k(:,:,i)));
        gamma_k(:,:,i) = B_k(:,:,i).*(C.');
    end
    
    %% 参数计算完毕，开始执行基于罚函数的连续凸近似方法
    reflect_vec = diag(reflect_mat);
    lamdba = 0.1;
    %控制迭代次数不高于100次
    for cnt_iter = 1:100
        reflect_vec_tmp = reflect_vec;%将上一次结果作为初始点
        lambda = lambda * 2;
        
        % 当每个IRS反射系数的模接近于1的时候，跳出循环
        if(sum(abs(abs(reflect_vec)-ones(scene.m_IRS,1))) < 0.01)
            break;
        end
    end
    
    % 将求得的近似反射向量的模置为1，得到结果
    reflect_mat = diag(exp(1j*angle(reflect_vec)));
end