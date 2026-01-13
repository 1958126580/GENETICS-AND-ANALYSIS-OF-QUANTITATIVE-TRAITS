# =============================================================================
# BayesAlphabet.jl - 贝叶斯字母表方法
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 21
# =============================================================================

"""
    bayesian_abc_models(y::AbstractVector, Z::AbstractMatrix;
                       prior_type::Symbol=:BayesA,
                       n_iter::Int=1000, burn_in::Int=100) -> PredictionResult

贝叶斯字母表模型 (BayesA/B/Cπ)。

# 方法
- :BayesA: 每个标记独立的t分布先验
- :BayesB: 混合先验（某些效应为0）
- :BayesCpi: 类似BayesB但π未知
"""
function bayesian_abc_models(y::AbstractVector, Z::AbstractMatrix;
                             prior_type::Symbol=:BayesA,
                             n_iter::Int=1000, burn_in::Int=100)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    Z_valid = Float64.(Z[valid_idx, :])
    
    n, m = size(Z_valid)
    
    # 初始化
    alpha = zeros(m)
    sigma2_e = var(y_valid) * 0.5
    sigma2_alpha = var(y_valid) * 0.5 / m
    
    # 存储
    alpha_samples = zeros(n_iter, m)
    
    for iter in 1:n_iter
        # Gibbs采样每个标记效应
        for j in 1:m
            z_j = Z_valid[:, j]
            resid = y_valid .- Z_valid * alpha .+ z_j .* alpha[j]
            
            # 条件后验
            ztz = dot(z_j, z_j)
            ztr = dot(z_j, resid)
            
            if prior_type == :BayesA
                var_post = 1 / (ztz/sigma2_e + 1/sigma2_alpha)
                mean_post = var_post * ztr / sigma2_e
                alpha[j] = mean_post + sqrt(var_post) * randn()
            elseif prior_type == :BayesB || prior_type == :BayesCpi
                pi = 0.9  # 简化
                if rand() < pi
                    alpha[j] = 0.0
                else
                    var_post = 1 / (ztz/sigma2_e + 1/sigma2_alpha)
                    mean_post = var_post * ztr / sigma2_e
                    alpha[j] = mean_post + sqrt(var_post) * randn()
                end
            end
        end
        
        # 更新sigma2_e
        resid = y_valid .- Z_valid * alpha
        sse = sum(resid.^2)
        sigma2_e = sse / rand(Chisq(n))
        
        alpha_samples[iter, :] = alpha
    end
    
    # 后验均值
    alpha_post = vec(mean(alpha_samples[(burn_in+1):end, :], dims=1))
    gebv = Z_valid * alpha_post
    
    accuracy = cor(gebv, y_valid)
    
    return PredictionResult(gebv, ["ind_$i" for i in valid_idx];
                           accuracy=accuracy, method=prior_type,
                           marker_effects=alpha_post)
end
