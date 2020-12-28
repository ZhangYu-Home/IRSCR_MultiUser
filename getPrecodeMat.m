%% 计算预编码矩阵
function precode_mat = getPrecodeMat(scene,g_AP_PUs,g_AP_SUs,decode_mat,weight_mat)
    %% 计算相关参数
    X_0 = zeros(scene.n_ante_AP,scene.n_ante_AP);
    Y = zeros(scene.n_data,scene.n_ante_AP,scene.n_SU);
    for i = 1:scene.n_SU
        X_0 = X_0 + g_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
        Y(:,:,i) = weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
    end
    [U,V] = eig(X_0);
    X_sqrt_0 = U*sqrt(V)*inv(U);
    
    X_k = zeros(scene.n_ante_AP,scene.n_ante_AP,scene.n_PU);
    X_sqrt_k = zeros(scene.n_ante_AP,scene.n_ante_AP,scene.n_PU);
    for i = 1:scene.n_PU
        X_k(:,:,i) = g_AP_PUs(:,:,i)'*g_AP_PUs(:,:,i);
        [U,V] = eig(X_k(:,:,i));
        X_sqrt_k(:,:,i) = U*sqrt(V)*inv(U);
    end
    
    
    
    cvx_begin quiet     %CVX solves convex problem
        variable F(scene.n_ante_AP,scene.n_data,scene.n_SU) complex
        objFunc = 0;
        consPower = 0;
        expression consInter(scene.n_PU,1);
        
        for i = 1:scene.n_SU
            objFunc = objFunc + square(norm(X_sqrt_0*F(:,:,i),'fro')) - trace(Y(:,:,i)*F(:,:,i)) - trace(Y(:,:,i)'*F(:,:,i)');
            consPower = consPower + square(norm(F(:,:,i),'fro'));
            for j = 1:scene.n_PU
                consInter(j) = consInter(j) + square(norm(X_sqrt_k(:,:,i)*F(:,:,i),'fro'));
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