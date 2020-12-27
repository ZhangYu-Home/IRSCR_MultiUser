%% 计算预编码矩阵
function precode_mat = getPrecodeMat(scene,g_AP_PUs,g_AP_SUs,decode_mat,weight_mat)
    %% 计算相关参数
    X_0 = zeros(scene.n_ante_AP,scene.n_ante_AP);
    Y = zeros(scene.n_data,scene.n_ante_AP,scene.n_SU);
    for i = 1:scene.n_SU
        X_0 = X_0 + g_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
        Y(:,:,i) = weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
    end
    
    X_k = zeros(scene.n_ante_AP,scene.n_ante_AP,scene.n_PU);
    for i = 1:scene.n_PU
        X_k(:,:,i) = g_AP_PUs(:,:,i)'*g_AP_PUs(:,:,i);
    end
    
    cvx_begin
    
    
    precode_mat = zeros(scene.n_ante_AP,scene.n_data,scene.n_SU);
    
end