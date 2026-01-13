# =============================================================================
# Association.jl - GWAS关联分析
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 20
# =============================================================================

"""
    gwas_single_marker_scan(y::AbstractVector, genotypes::AbstractMatrix;
                           covariates::Union{Matrix, Nothing}=nothing) -> GWASResult

单标记GWAS扫描。
"""
function gwas_single_marker_scan(y::AbstractVector, genotypes::AbstractMatrix;
                                 covariates::Union{AbstractMatrix, Nothing}=nothing)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    geno_valid = Float64.(genotypes[valid_idx, :])
    n = length(y_valid)
    n_snps = size(geno_valid, 2)
    
    # 协变量
    X0 = covariates === nothing ? ones(n, 1) : hcat(ones(n), covariates[valid_idx, :])
    
    pvalues = zeros(n_snps)
    betas = zeros(n_snps)
    ses = zeros(n_snps)
    mafs = zeros(n_snps)
    
    for j in 1:n_snps
        g = geno_valid[:, j]
        
        # MAF
        p = mean(g) / 2
        mafs[j] = min(p, 1-p)
        
        # 回归
        X = hcat(X0, g)
        try
            beta = X \ y_valid
            y_hat = X * beta
            res = y_valid .- y_hat
            MSE = sum(res.^2) / (n - size(X, 2))
            
            XtX_inv = inv(X' * X)
            se = sqrt(MSE * XtX_inv[end, end])
            t_stat = beta[end] / se
            pval = 2 * ccdf(TDist(n - size(X, 2)), abs(t_stat))
            
            betas[j] = beta[end]
            ses[j] = se
            pvalues[j] = pval
        catch
            pvalues[j] = 1.0
        end
    end
    
    # 基因组控制
    lambda_gc = genomic_control_lambda(pvalues)
    
    return GWASResult(
        ["SNP_$j" for j in 1:n_snps],
        fill("Chr1", n_snps),
        Float64.(1:n_snps),
        pvalues, betas, ses, mafs;
        lambda_gc=lambda_gc, method=:GLM
    )
end

"""
    gwas_mlm_fit(y::AbstractVector, genotypes::AbstractMatrix, 
                K::AbstractMatrix; PC::Union{Matrix, Nothing}=nothing) -> GWASResult

混合线性模型GWAS (EMMA/EMMAX风格)。
"""
function gwas_mlm_fit(y::AbstractVector, genotypes::AbstractMatrix,
                     K::AbstractMatrix; PC::Union{AbstractMatrix, Nothing}=nothing)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    geno_valid = Float64.(genotypes[valid_idx, :])
    K_valid = K[valid_idx, valid_idx]
    n = length(y_valid)
    n_snps = size(geno_valid, 2)
    
    # 估计方差组分
    X0 = PC === nothing ? ones(n, 1) : hcat(ones(n), PC[valid_idx, :])
    vc = reml_ai_algorithm(y_valid, X0, K_valid)
    
    delta = vc.sigma2_e / vc.sigma2_a
    H = K_valid + delta * I(n)
    H_inv = inv(H)
    
    pvalues = zeros(n_snps)
    betas = zeros(n_snps)
    ses = zeros(n_snps)
    mafs = zeros(n_snps)
    
    for j in 1:n_snps
        g = geno_valid[:, j]
        mafs[j] = min(mean(g)/2, 1 - mean(g)/2)
        
        X = hcat(X0, g)
        
        # GLS
        XtHinv = X' * H_inv
        XtHinvX = XtHinv * X
        XtHinvy = XtHinv * y_valid
        
        try
            beta = XtHinvX \ XtHinvy
            y_hat = X * beta
            res = y_valid .- y_hat
            
            se = sqrt(diag(inv(XtHinvX))[end] * vc.sigma2_e)
            t_stat = beta[end] / se
            pval = 2 * ccdf(TDist(n - size(X, 2)), abs(t_stat))
            
            betas[j] = beta[end]
            ses[j] = se
            pvalues[j] = pval
        catch
            pvalues[j] = 1.0
        end
    end
    
    lambda_gc = genomic_control_lambda(pvalues)
    
    return GWASResult(
        ["SNP_$j" for j in 1:n_snps],
        fill("Chr1", n_snps),
        Float64.(1:n_snps),
        pvalues, betas, ses, mafs;
        lambda_gc=lambda_gc, method=:MLM
    )
end

"""
    genomic_control_lambda(pvalues::Vector{Float64}) -> Float64

计算基因组控制膨胀因子λ_GC。
"""
function genomic_control_lambda(pvalues::Vector{Float64})
    valid_p = filter(p -> 0 < p < 1, pvalues)
    isempty(valid_p) && return 1.0
    
    chi2_stats = [quantile(Chisq(1), 1 - p) for p in valid_p]
    median_chi2 = median(chi2_stats)
    expected_median = 0.4549  # median of χ²(1)
    
    return median_chi2 / expected_median
end

"""
    structured_association_mapping(y::AbstractVector, genotypes::AbstractMatrix,
                                  Q::AbstractMatrix; K::Union{AbstractMatrix,Nothing}=nothing) -> GWASResult

结构化关联分析 (Q+K模型)。

# 模型
y = Xβ + Qv + Zu + ε

其中:
- Q: 群体结构矩阵（如STRUCTURE或PCA前几个PC）
- K: 亲缘关系矩阵（可选）

# 参数
- `y`: 表型向量
- `genotypes`: 基因型矩阵
- `Q`: 群体结构矩阵
- `K`: 亲缘关系矩阵（可选，如果提供则使用Q+K模型）

# 中文说明
用于存在群体分层的关联分析，通过纳入群体结构和亲缘关系矩阵控制假阳性。

# 公式来源
Yu et al. (2006) Nature Genetics, Walsh 2nd Ed Chapter 20
"""
function structured_association_mapping(y::AbstractVector, genotypes::AbstractMatrix,
                                        Q::AbstractMatrix; 
                                        K::Union{AbstractMatrix,Nothing}=nothing)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    geno_valid = Float64.(genotypes[valid_idx, :])
    Q_valid = Float64.(Q[valid_idx, :])
    n = length(y_valid)
    n_snps = size(geno_valid, 2)
    
    # 构建协变量矩阵（截距 + Q矩阵）
    X0 = hcat(ones(n), Q_valid)
    
    pvalues = zeros(n_snps)
    betas = zeros(n_snps)
    ses = zeros(n_snps)
    mafs = zeros(n_snps)
    
    if K !== nothing
        # Q+K模型：使用混合模型
        K_valid = K[valid_idx, valid_idx]
        
        # 估计方差组分
        vc = reml_ai_algorithm(y_valid, X0, K_valid; max_iter=50)
        
        if vc.sigma2_a > 0
            delta = vc.sigma2_e / vc.sigma2_a
            H = K_valid + delta * I(n)
        else
            H = I(n)
        end
        
        H_inv = inv(H + 1e-6 * I(n))
        
        for j in 1:n_snps
            g = geno_valid[:, j]
            mafs[j] = min(mean(g)/2, 1 - mean(g)/2)
            
            X = hcat(X0, g)
            
            try
                XtHinv = X' * H_inv
                XtHinvX = XtHinv * X
                XtHinvy = XtHinv * y_valid
                
                beta = XtHinvX \ XtHinvy
                
                se = sqrt(abs(diag(inv(XtHinvX))[end]) * vc.sigma2_e)
                t_stat = beta[end] / max(se, 1e-10)
                pval = 2 * ccdf(TDist(max(n - size(X, 2), 1)), abs(t_stat))
                
                betas[j] = beta[end]
                ses[j] = se
                pvalues[j] = pval
            catch
                pvalues[j] = 1.0
            end
        end
    else
        # 仅Q模型：普通线性回归
        for j in 1:n_snps
            g = geno_valid[:, j]
            mafs[j] = min(mean(g)/2, 1 - mean(g)/2)
            
            X = hcat(X0, g)
            
            try
                beta = X \ y_valid
                y_hat = X * beta
                res = y_valid .- y_hat
                MSE = sum(res.^2) / (n - size(X, 2))
                
                XtX_inv = inv(X' * X)
                se = sqrt(MSE * XtX_inv[end, end])
                t_stat = beta[end] / se
                pval = 2 * ccdf(TDist(n - size(X, 2)), abs(t_stat))
                
                betas[j] = beta[end]
                ses[j] = se
                pvalues[j] = pval
            catch
                pvalues[j] = 1.0
            end
        end
    end
    
    lambda_gc = genomic_control_lambda(pvalues)
    
    return GWASResult(
        ["SNP_$j" for j in 1:n_snps],
        fill("Chr1", n_snps),
        Float64.(1:n_snps),
        pvalues, betas, ses, mafs;
        lambda_gc=lambda_gc, method=K !== nothing ? :QK : :Q
    )
end

"""
    bonferroni_threshold(n_tests::Int; alpha::Float64=0.05) -> Float64

Bonferroni校正阈值。
"""
bonferroni_threshold(n_tests::Int; alpha::Float64=0.05) = alpha / n_tests

"""
    fdr_correction(pvalues::Vector{Float64}; alpha::Float64=0.05) -> NamedTuple

Benjamini-Hochberg FDR校正。

# 返回值
- `significant`: 显著的SNP索引
- `q_values`: 调整后的p值（q值）
- `threshold`: FDR阈值
"""
function fdr_correction(pvalues::Vector{Float64}; alpha::Float64=0.05)
    n = length(pvalues)
    sorted_idx = sortperm(pvalues)
    sorted_p = pvalues[sorted_idx]
    
    # BH阈值
    thresholds = [(i / n) * alpha for i in 1:n]
    
    # 找到最大的k使得p(k) <= k/n * α
    significant_k = 0
    for k in 1:n
        if sorted_p[k] <= thresholds[k]
            significant_k = k
        end
    end
    
    # q值计算
    q_values = zeros(n)
    q_values[sorted_idx[n]] = sorted_p[n]
    for i in (n-1):-1:1
        q_values[sorted_idx[i]] = min(sorted_p[i] * n / i, q_values[sorted_idx[i+1]])
    end
    
    significant = significant_k > 0 ? sorted_idx[1:significant_k] : Int[]
    
    return (significant=significant, q_values=q_values, 
            threshold=significant_k > 0 ? thresholds[significant_k] : 0.0)
end

"""
    ld_clumping(pvalues::Vector{Float64}, ld_matrix::AbstractMatrix;
               p_threshold::Float64=5e-8, r2_threshold::Float64=0.2,
               window_kb::Float64=500.0) -> Vector{Int}

LD clumping选择独立显著SNP。

# 算法
1. 选择最显著的SNP作为lead SNP
2. 移除与lead SNP LD > r²阈值的SNP
3. 重复直到没有显著SNP
"""
function ld_clumping(pvalues::Vector{Float64}, ld_matrix::AbstractMatrix;
                     p_threshold::Float64=5e-8, r2_threshold::Float64=0.2,
                     window_kb::Float64=500.0)
    n = length(pvalues)
    
    # 初始化
    available = trues(n)
    lead_snps = Int[]
    
    while true
        # 找到最显著的可用SNP
        best_idx = 0
        best_p = Inf
        
        for i in 1:n
            if available[i] && pvalues[i] < p_threshold && pvalues[i] < best_p
                best_p = pvalues[i]
                best_idx = i
            end
        end
        
        if best_idx == 0
            break
        end
        
        push!(lead_snps, best_idx)
        available[best_idx] = false
        
        # 移除与lead SNP高LD的SNP
        for j in 1:n
            if available[j]
                r2 = ld_matrix[best_idx, j]^2
                if r2 > r2_threshold
                    available[j] = false
                end
            end
        end
    end
    
    return lead_snps
end

