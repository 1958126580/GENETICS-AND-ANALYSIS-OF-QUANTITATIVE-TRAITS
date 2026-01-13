# =============================================================================
# Moments.jl - 矩与累积量分析
# =============================================================================
#
# 实现分布的矩、累积量和正态性变换
#
# 对应章节: Walsh 2nd Ed, Chapter 2 (Properties of Distributions)
# =============================================================================

"""
    moments_around_mean(data::AbstractVector, k::Int) -> Float64

计算数据的第k阶中心矩。

# 数学定义
```math
\\mu_k = E[(X - \\mu)^k] = \\frac{1}{n}\\sum_{i=1}^{n}(x_i - \\bar{x})^k
```

# 参数
- `data::AbstractVector`: 数据向量
- `k::Int`: 矩的阶数（k ≥ 1）

# 返回值
- `Float64`: 第k阶中心矩

# 中文说明
中心矩描述数据围绕均值的分布特征：
- μ₁ = 0（恒为零）
- μ₂ = 方差（σ²）
- μ₃ 与偏度相关
- μ₄ 与峰度相关

# 公式来源
Walsh 2nd Ed, Eq 2.3-2.5

# 示例
```julia
data = randn(1000)
moments_around_mean(data, 2)  # ≈ 1.0 (方差)
moments_around_mean(data, 3)  # ≈ 0.0 (正态分布偏度)
moments_around_mean(data, 4)  # ≈ 3.0 (正态分布峰度)
```

参见: [`cumulants_from_moments`](@ref), [`skewness`](@ref), [`kurtosis`](@ref)
"""
function moments_around_mean(data::AbstractVector, k::Int)
    if k < 1
        throw(ArgumentError("矩阶数必须 ≥ 1"))
    end
    
    # 移除缺失值
    valid = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)), data)
    n = length(valid)
    
    if n == 0
        throw(ArgumentError("没有有效的观测值"))
    end
    
    x = Float64.(valid)
    mu = mean(x)
    
    return mean((x .- mu).^k)
end

"""
    raw_moments(data::AbstractVector, k::Int) -> Float64

计算数据的第k阶原始矩（关于原点）。

# 数学定义
```math
m_k = E[X^k] = \\frac{1}{n}\\sum_{i=1}^{n}x_i^k
```

# 中文说明
原始矩是关于原点（0）计算的矩：
- m₁ = 均值
- m₂ = E[X²]
"""
function raw_moments(data::AbstractVector, k::Int)
    if k < 1
        throw(ArgumentError("矩阶数必须 ≥ 1"))
    end
    
    valid = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)), data)
    x = Float64.(valid)
    
    return mean(x.^k)
end

"""
    skewness(data::AbstractVector) -> Float64

计算数据的偏度（标准化三阶矩）。

# 数学定义
```math
\\gamma_1 = \\frac{\\mu_3}{\\sigma^3} = \\frac{E[(X-\\mu)^3]}{(E[(X-\\mu)^2])^{3/2}}
```

# 中文说明
偏度衡量分布的不对称程度：
- γ₁ = 0: 对称分布（如正态分布）
- γ₁ > 0: 右偏（长尾在右）
- γ₁ < 0: 左偏（长尾在左）

# 示例
```julia
skewness(randn(1000))  # ≈ 0.0
skewness(rand(1000).^2)  # > 0 (右偏)
```
"""
function skewness(data::AbstractVector)
    mu_2 = moments_around_mean(data, 2)
    mu_3 = moments_around_mean(data, 3)
    
    if mu_2 ≈ 0
        return NaN
    end
    
    return mu_3 / (mu_2^1.5)
end

"""
    kurtosis(data::AbstractVector; excess::Bool=true) -> Float64

计算数据的峰度。

# 数学定义
```math
\\gamma_2 = \\frac{\\mu_4}{\\sigma^4} - 3 \\quad \\text{(超额峰度)}
```

# 参数
- `data`: 数据向量
- `excess`: 是否返回超额峰度（默认true）

# 中文说明
峰度衡量分布尾部的厚度：
- γ₂ = 0: 正态分布（mesokurtic）
- γ₂ > 0: 厚尾分布（leptokurtic）
- γ₂ < 0: 薄尾分布（platykurtic）

注意：正态分布的原始峰度 = 3，超额峰度 = 0

# 示例
```julia
kurtosis(randn(10000))  # ≈ 0.0 (正态)
kurtosis(rand(10000))   # ≈ -1.2 (均匀分布)
```
"""
function kurtosis(data::AbstractVector; excess::Bool=true)
    mu_2 = moments_around_mean(data, 2)
    mu_4 = moments_around_mean(data, 4)
    
    if mu_2 ≈ 0
        return NaN
    end
    
    raw_kurt = mu_4 / (mu_2^2)
    
    return excess ? raw_kurt - 3.0 : raw_kurt
end

"""
    cumulants_from_moments(moments::Vector{Float64}) -> Vector{Float64}

从原始矩计算累积量。

# 数学关系
```math
\\begin{aligned}
\\kappa_1 &= \\mu'_1 = \\mu \\quad \\text{(均值)} \\\\
\\kappa_2 &= \\mu_2 = \\sigma^2 \\quad \\text{(方差)} \\\\
\\kappa_3 &= \\mu_3 \\quad \\text{(三阶中心矩)} \\\\
\\kappa_4 &= \\mu_4 - 3\\mu_2^2 \\quad \\text{(与峰度相关)}
\\end{aligned}
```

# 参数
- `moments`: 前k阶原始矩的向量 [m₁, m₂, m₃, m₄]

# 返回值
- `Vector{Float64}`: 对应的累积量 [κ₁, κ₂, κ₃, κ₄]

# 中文说明
累积量（Cumulant）是矩的替代表示，具有重要的统计性质：
- 独立随机变量之和的累积量等于各累积量之和
- 正态分布的所有 κₖ (k ≥ 3) 都等于0

# 公式来源
Walsh 2nd Ed, Eq 2.8-2.11

# 示例
```julia
# 从数据计算累积量
data = randn(1000)
m = [raw_moments(data, k) for k in 1:4]
kappa = cumulants_from_moments(m)
```
"""
function cumulants_from_moments(moments::Vector{Float64})
    n = length(moments)
    if n < 1
        throw(ArgumentError("至少需要一阶矩"))
    end
    
    kappa = zeros(n)
    
    # κ₁ = μ'₁ (均值)
    kappa[1] = moments[1]
    
    if n >= 2
        # κ₂ = μ'₂ - μ'₁² (方差)
        kappa[2] = moments[2] - moments[1]^2
    end
    
    if n >= 3
        # κ₃ = μ'₃ - 3μ'₁μ'₂ + 2μ'₁³
        kappa[3] = moments[3] - 3*moments[1]*moments[2] + 2*moments[1]^3
    end
    
    if n >= 4
        # κ₄ = μ'₄ - 4μ'₁μ'₃ - 3μ'₂² + 12μ'₁²μ'₂ - 6μ'₁⁴
        kappa[4] = moments[4] - 4*moments[1]*moments[3] - 3*moments[2]^2 + 
                   12*moments[1]^2*moments[2] - 6*moments[1]^4
    end
    
    return kappa
end

"""
    multivariate_normal_dens(X::AbstractVector, mu::AbstractVector, 
                            Sigma::AbstractMatrix; log::Bool=true) -> Float64

计算多元正态分布的（对数）密度。

# 数学定义
```math
f(\\mathbf{x}) = \\frac{1}{(2\\pi)^{p/2}|\\Sigma|^{1/2}} 
\\exp\\left(-\\frac{1}{2}(\\mathbf{x}-\\boldsymbol{\\mu})^T \\Sigma^{-1} (\\mathbf{x}-\\boldsymbol{\\mu})\\right)
```

# 参数
- `X`: 观测向量（p维）
- `mu`: 均值向量（p维）
- `Sigma`: 协方差矩阵（p×p，必须正定）
- `log`: 是否返回对数密度（默认true）

# 中文说明
多元正态分布是数量遗传学的核心分布，用于：
- 多性状分析
- 混合模型似然计算
- 贝叶斯推断

使用对数形式可避免数值下溢。

# 公式来源
Walsh 2nd Ed, Eq 2.15

# 示例
```julia
mu = [0.0, 0.0]
Sigma = [1.0 0.5; 0.5 1.0]
x = [0.5, 0.3]
log_density = multivariate_normal_dens(x, mu, Sigma)
density = exp(log_density)
```
"""
function multivariate_normal_dens(X::AbstractVector, mu::AbstractVector,
                                  Sigma::AbstractMatrix; log::Bool=true)
    p = length(X)
    
    if length(mu) != p
        throw(ArgumentError("X和mu维度不匹配"))
    end
    if size(Sigma) != (p, p)
        throw(ArgumentError("Sigma维度与X不匹配"))
    end
    
    # Cholesky分解（更稳定）
    C = cholesky(Symmetric(Sigma))
    
    # 计算 (x - mu)' Σ^(-1) (x - mu)
    diff = X .- mu
    z = C.L \ diff  # 求解 L*z = diff
    mahalanobis = dot(z, z)
    
    # log|Σ| = 2 * sum(log(diag(L)))
    log_det = 2 * sum(log.(diag(C.L)))
    
    # 对数密度
    log_dens = -0.5 * (p * log(2π) + log_det + mahalanobis)
    
    return log ? log_dens : exp(log_dens)
end

"""
    mixture_model_likelihood(y::AbstractVector, 
                            components::Vector{<:Distribution};
                            weights::AbstractVector=ones(length(components))/length(components)) -> Float64

计算混合模型的（对数）似然。

# 数学定义
```math
L(\\mathbf{y}) = \\prod_{i=1}^{n} \\left[ \\sum_{j=1}^{K} \\pi_j f_j(y_i) \\right]
```

# 参数
- `y`: 观测值向量
- `components`: 分布成分向量（如正态分布列表）
- `weights`: 混合权重（默认等权重，和为1）

# 中文说明
混合模型假设观测来自多个亚群体，每个亚群体有不同的分布。
在QTL分析中，不同基因型对应不同的正态分布成分。

# 公式来源
Walsh 2nd Ed, Chapter 18

# 示例
```julia
# 双成分正态混合
comp1 = Normal(0.0, 1.0)
comp2 = Normal(3.0, 1.0)
y = [randn() for _ in 1:50]
append!(y, [randn() + 3.0 for _ in 1:50])

ll = mixture_model_likelihood(y, [comp1, comp2], weights=[0.5, 0.5])
```
"""
function mixture_model_likelihood(y::AbstractVector,
                                  components::Vector{<:Distribution};
                                  weights::AbstractVector=ones(length(components))/length(components))
    K = length(components)
    if length(weights) != K
        throw(ArgumentError("权重数量与成分数量不匹配"))
    end
    
    # 归一化权重
    w = weights ./ sum(weights)
    
    # 过滤有效观测
    valid_y = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)), y)
    n = length(valid_y)
    
    if n == 0
        return -Inf
    end
    
    # 计算对数似然
    log_lik = 0.0
    for yi in valid_y
        # 计算 sum(w_k * f_k(y))
        log_densities = [log(w[k]) + logpdf(components[k], yi) for k in 1:K]
        log_lik += logsumexp(log_densities)
    end
    
    return log_lik
end

"""
    transformation_normal(y::AbstractVector; 
                         method::Symbol=:boxcox,
                         lambda::Union{Float64, Nothing}=nothing) -> NamedTuple

对数据进行正态化变换。

# 方法
- `:boxcox`: Box-Cox变换（最常用）
- `:log`: 对数变换（λ=0的特例）
- `:sqrt`: 平方根变换（λ=0.5的特例）
- `:rank`: 秩逆正态变换

# Box-Cox变换
```math
y^{(\\lambda)} = \\begin{cases}
\\frac{y^\\lambda - 1}{\\lambda} & \\lambda \\neq 0 \\\\
\\log(y) & \\lambda = 0
\\end{cases}
```

# 中文说明
当表型数据不满足正态性假设时，需要进行变换。
Box-Cox变换通过最大似然估计选择最优的λ参数。

常见变换：
- λ = 1: 无变换
- λ = 0: 对数变换
- λ = 0.5: 平方根变换
- λ = -1: 倒数变换

# 公式来源
Walsh 2nd Ed, Chapter 14 (Scale)

# 示例
```julia
# 自动Box-Cox变换
y_raw = rand(100).^3  # 右偏数据
result = transformation_normal(y_raw, method=:boxcox)
result.y_transformed  # 变换后的数据
result.lambda         # 最优lambda
```
"""
function transformation_normal(y::AbstractVector;
                               method::Symbol=:boxcox,
                               lambda::Union{Float64, Nothing}=nothing)
    valid_y = filter(v -> !ismissing(v) && !(v isa Number && isnan(v)) && v > 0, y)
    
    if isempty(valid_y)
        throw(ArgumentError("没有有效的正值观测"))
    end
    
    y_float = Float64.(valid_y)
    
    if method == :log
        y_transformed = log.(y_float)
        optimal_lambda = 0.0
        
    elseif method == :sqrt
        y_transformed = sqrt.(y_float)
        optimal_lambda = 0.5
        
    elseif method == :rank
        # 秩逆正态变换
        n = length(y_float)
        ranks = ordinalrank(y_float)
        # 逆正态变换
        y_transformed = [quantile(Normal(), (r - 0.5) / n) for r in ranks]
        optimal_lambda = NaN
        
    elseif method == :boxcox
        if lambda === nothing
            # 通过最大似然估计最优lambda
            optimal_lambda = optimize_boxcox_lambda(y_float)
        else
            optimal_lambda = lambda
        end
        y_transformed = boxcox_transform(y_float, optimal_lambda)
        
    else
        throw(ArgumentError("未知的变换方法: $method"))
    end
    
    # 计算变换后的正态性检验
    sw_stat = shapiro_wilk_approx(y_transformed)
    
    return (
        y_transformed = y_transformed,
        y_original = y_float,
        lambda = optimal_lambda,
        method = method,
        normality_stat = sw_stat
    )
end

# Box-Cox变换
function boxcox_transform(y::AbstractVector, lambda::Float64)
    if abs(lambda) < 1e-10
        return log.(y)
    else
        return (y.^lambda .- 1) ./ lambda
    end
end

# 优化Box-Cox的lambda参数
function optimize_boxcox_lambda(y::AbstractVector)
    # 简单网格搜索
    lambdas = range(-2, 2, length=41)
    best_lambda = 1.0
    best_loglik = -Inf
    
    n = length(y)
    
    for lam in lambdas
        y_t = boxcox_transform(y, lam)
        
        # 正态分布对数似然（忽略常数项）
        sigma2 = var(y_t)
        if sigma2 > 0
            loglik = -n/2 * log(sigma2) + (lam - 1) * sum(log.(y))
            if loglik > best_loglik
                best_loglik = loglik
                best_lambda = lam
            end
        end
    end
    
    return best_lambda
end

# Shapiro-Wilk近似统计量
function shapiro_wilk_approx(y::AbstractVector)
    n = length(y)
    if n < 3
        return NaN
    end
    
    y_sorted = sort(y)
    y_mean = mean(y)
    
    # 简化的W统计量近似
    ss = sum((y .- y_mean).^2)
    
    if ss ≈ 0
        return 1.0
    end
    
    # 使用正态分位数近似
    m = [quantile(Normal(), (i - 0.375) / (n + 0.25)) for i in 1:n]
    a = m ./ norm(m)
    
    b = dot(a, y_sorted)
    W = b^2 / ss
    
    return W
end

# 秩函数（如果未导入）
function ordinalrank(x::AbstractVector)
    n = length(x)
    perm = sortperm(x)
    ranks = zeros(Int, n)
    for (rank, idx) in enumerate(perm)
        ranks[idx] = rank
    end
    return ranks
end
