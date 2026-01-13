# =============================================================================
# REML.jl - 限制性最大似然方差组分估计
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 10
# =============================================================================

"""
    reml_ai_algorithm(y::AbstractVector, X::AbstractMatrix, K::AbstractMatrix;
                      max_iter::Int=100, tol::Float64=1e-6) -> VarianceComponents

AI-REML算法估计方差组分。

# 模型
y = Xβ + u + ε
Var(u) = Kσ²_a, Var(ε) = Iσ²_e

# 中文说明
Average Information REML比EM-REML收敛更快。
"""
function reml_ai_algorithm(y::AbstractVector, X::AbstractMatrix, K::AbstractMatrix;
                           max_iter::Int=100, tol::Float64=1e-6)
    n = length(y)
    p = size(X, 2)
    
    # 初始值
    var_y = var(y)
    sigma2_a = var_y * 0.5
    sigma2_e = var_y * 0.5
    
    converged = false
    iter = 0
    
    for i in 1:max_iter
        iter = i
        
        # V = Kσ²_a + Iσ²_e
        V = K * sigma2_a + I(n) * sigma2_e
        
        # 尝试Cholesky分解
        V_chol = try
            cholesky(Symmetric(V))
        catch
            V = V + 1e-6 * I(n)
            cholesky(Symmetric(V))
        end
        
        V_inv = inv(V_chol)
        
        # P = V⁻¹ - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹
        XtVinv = X' * V_inv
        XtVinvX_inv = inv(XtVinv * X)
        P = V_inv - V_inv * X * XtVinvX_inv * XtVinv
        
        # REML似然的一阶导数
        Py = P * y
        
        dL_da = -0.5 * (tr(P * K) - Py' * K * Py)
        dL_de = -0.5 * (tr(P) - Py' * Py)
        
        # AI矩阵（二阶导数近似）
        PKPy = P * K * Py
        PPy = P * Py
        
        AI_aa = 0.5 * Py' * K * P * K * Py
        AI_ae = 0.5 * Py' * K * P * Py
        AI_ee = 0.5 * Py' * P * Py
        
        AI = [AI_aa AI_ae; AI_ae AI_ee]
        grad = [dL_da; dL_de]
        
        # Newton-Raphson更新
        delta = AI \ grad
        
        sigma2_a_new = sigma2_a + delta[1]
        sigma2_e_new = sigma2_e + delta[2]
        
        # 确保方差为正
        sigma2_a_new = max(sigma2_a_new, 1e-6)
        sigma2_e_new = max(sigma2_e_new, 1e-6)
        
        # 检查收敛
        if abs(sigma2_a_new - sigma2_a) < tol && abs(sigma2_e_new - sigma2_e) < tol
            converged = true
            sigma2_a = sigma2_a_new
            sigma2_e = sigma2_e_new
            break
        end
        
        sigma2_a = sigma2_a_new
        sigma2_e = sigma2_e_new
    end
    
    return VarianceComponents(sigma2_a, sigma2_e; 
                             method=:REML, converged=converged, n_iterations=iter)
end

"""
    gibbs_sampler_vc(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix,
                    A::AbstractMatrix; n_iter::Int=1000, burn_in::Int=100) -> NamedTuple

Gibbs采样估计方差组分。

# 中文说明
贝叶斯MCMC方法，适用于复杂模型和小样本。
"""
function gibbs_sampler_vc(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix,
                          A::AbstractMatrix; n_iter::Int=1000, burn_in::Int=100)
    n = length(y)
    p = size(X, 2)
    q = size(Z, 2)
    
    # 初始化
    beta = X \ y
    u = zeros(q)
    sigma2_a = var(y) * 0.5
    sigma2_e = var(y) * 0.5
    
    # 存储样本
    samples_a = zeros(n_iter)
    samples_e = zeros(n_iter)
    
    A_inv = inv(A + 1e-6 * I(q))
    
    for iter in 1:n_iter
        # 采样β
        resid = y - Z * u
        XtX = X' * X
        beta = XtX \ (X' * resid) + randn(p) .* sqrt.(diag(inv(XtX)) .* sigma2_e)
        
        # 采样u
        resid_u = y - X * beta
        lambda = sigma2_e / sigma2_a
        lhs = Z' * Z + lambda * A_inv
        rhs = Z' * resid_u
        u_mean = lhs \ rhs
        u = u_mean + randn(q) .* sqrt.(diag(inv(lhs)) .* sigma2_e)
        
        # 采样σ²_e
        resid_full = y - X * beta - Z * u
        sse = sum(resid_full.^2)
        sigma2_e = sse / rand(Chisq(n - p))
        
        # 采样σ²_a
        ssu = u' * A_inv * u
        sigma2_a = ssu / rand(Chisq(q))
        
        samples_a[iter] = sigma2_a
        samples_e[iter] = sigma2_e
    end
    
    # 后验均值（移除burn-in）
    post_sigma2_a = mean(samples_a[(burn_in+1):end])
    post_sigma2_e = mean(samples_e[(burn_in+1):end])
    
    return (sigma2_a=post_sigma2_a, sigma2_e=post_sigma2_e,
            samples_a=samples_a, samples_e=samples_e)
end
