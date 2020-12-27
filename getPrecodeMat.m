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
    
    objFunc = 0;
    consPower = 0;
    consInter = zeros(scene.n_PU,1);
    
    size(X_0)
    scene.n_ante_AP
    scene.n_data
    cvx_solver
    cvx_begin      %CVX solves convex problem
        variable F(scene.n_ante_AP,scene.n_data,scene.n_SU) complex
        for i = 1:scene.n_SU
            objFunc = objFunc + trace(F(:,:,i)'*X_0*F(:,:,i)) - trace(Y(:,:,i)*F(:,:,i)) - trace(Y(:,:,i)'*F(:,:,i)');
            consPower = consPower + trace(F(:,:,i)'*F(:,:,i));
            for j = 1:scene.n_PU
                consInter(j) = consInter(j) + trace(F(:,:,i)'*X_k(:,:,j)*F(:,:,i));
            end
        end
        minimize(objFunc)	%目标函数
        subject to
            %总功率约束
            consPower <= scene.max_pow;
            for j = 1:scene.n_PU
                consInter(j) <= scene.leak_pow;
            end
    cvx_end
    precode_mat = F;
end