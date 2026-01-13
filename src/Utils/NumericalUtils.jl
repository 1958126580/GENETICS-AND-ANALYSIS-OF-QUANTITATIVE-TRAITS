# =============================================================================
# NumericalUtils.jl - 数值计算工具函数
# =============================================================================
#
# 提供数值稳定的计算函数，防止溢出/下溢
#
# =============================================================================

"""
    safe_log(x::Real) -> Float64

安全的对数函数，防止log(0)产生-Inf。

# 参数
- `x::Real`: 输入值

# 返回值
- `Float64`: log(max(x, MIN_PROBABILITY))

# 中文说明
安全对数函数，当输入趋近于0时返回log(MIN_PROBABILITY)而非-Inf。
用于似然计算中防止数值下溢。

# 示例
```julia
safe_log(0.0)     # ≈ -690.8 而非 -Inf
safe_log(1e-400)  # ≈ -690.8 而非 -Inf
safe_log(0.5)     # ≈ -0.693
```
"""
function safe_log(x::Real)
    return log(max(x, MIN_PROBABILITY))
end

"""
    safe_exp(x::Real) -> Float64

安全的指数函数，防止exp溢出。

# 参数
- `x::Real`: 输入值

# 返回值
- `Float64`: 受限的exp(x)值

# 中文说明
安全指数函数，对输入值进行裁剪以防止溢出。
"""
function safe_exp(x::Real)
    # 裁剪输入以防止溢出
    x_clipped = clamp(x, -700.0, 700.0)
    return exp(x_clipped)
end

"""
    logit(p::Real) -> Float64

Logit变换: log(p/(1-p))

# 参数
- `p::Real`: 概率值，0 < p < 1

# 返回值
- `Float64`: logit(p)

# 中文说明
将概率值[0,1]映射到实数域(-∞, +∞)。
常用于参数变换以消除边界约束。

# 公式
logit(p) = log(p / (1 - p))

# 示例
```julia
logit(0.5)  # 0.0
logit(0.9)  # ≈ 2.197
```
"""
function logit(p::Real)
    p_safe = clamp(p, MIN_PROBABILITY, MAX_PROBABILITY)
    return log(p_safe / (1.0 - p_safe))
end

"""
    expit(x::Real) -> Float64

Logistic函数（expit），logit的逆变换。

# 公式
expit(x) = 1 / (1 + exp(-x)) = exp(x) / (1 + exp(x))

# 中文说明
将实数映射回概率[0,1]区间。
使用数值稳定的实现避免溢出。
"""
function expit(x::Real)
    if x >= 0
        z = exp(-x)
        return 1.0 / (1.0 + z)
    else
        z = exp(x)
        return z / (1.0 + z)
    end
end

"""
    logsumexp(x::AbstractVector{<:Real}) -> Float64

Log-Sum-Exp技巧，数值稳定地计算 log(sum(exp(x)))。

# 中文说明
用于在对数空间中将多个概率相加。
避免直接计算exp导致的溢出问题。

# 公式
logsumexp(x) = max(x) + log(sum(exp(x - max(x))))

# 示例
```julia
logsumexp([-1000, -1001, -1002])  # ≈ -999.59（而非-Inf）
```
"""
function logsumexp(x::AbstractVector{<:Real})
    if isempty(x)
        return -Inf
    end
    max_x = maximum(x)
    if isinf(max_x)
        return max_x
    end
    return max_x + log(sum(exp.(x .- max_x)))
end

"""
    logsumexp(a::Real, b::Real) -> Float64

两个值的Log-Sum-Exp。
"""
function logsumexp(a::Real, b::Real)
    if a > b
        return a + log1p(exp(b - a))
    else
        return b + log1p(exp(a - b))
    end
end

"""
    logdiffexp(a::Real, b::Real) -> Float64

计算 log(exp(a) - exp(b))，要求 a > b。

# 中文说明
用于在对数空间中计算概率差。
"""
function logdiffexp(a::Real, b::Real)
    if a <= b
        return -Inf
    end
    return a + log1p(-exp(b - a))
end

"""
    normalize_log_probabilities(log_probs::AbstractVector{<:Real}) -> Vector{Float64}

将对数概率归一化为和为1的常规概率。

# 中文说明
将对数概率转换为归一化的概率分布。
使用logsumexp保持数值稳定。

# 示例
```julia
log_probs = [-10.0, -11.0, -12.0]
probs = normalize_log_probabilities(log_probs)
sum(probs)  # ≈ 1.0
```
"""
function normalize_log_probabilities(log_probs::AbstractVector{<:Real})
    log_total = logsumexp(log_probs)
    return exp.(log_probs .- log_total)
end

"""
    is_positive_definite(A::AbstractMatrix) -> Bool

检查矩阵是否正定。

# 中文说明
通过尝试Cholesky分解来判断矩阵正定性。
"""
function is_positive_definite(A::AbstractMatrix)
    try
        cholesky(Symmetric(A))
        return true
    catch
        return false
    end
end

"""
    make_positive_definite(A::AbstractMatrix; epsilon::Float64=1e-6) -> Matrix{Float64}

将矩阵调整为正定矩阵。

# 方法
通过特征值分解，将负特征值替换为小正值。

# 中文说明
某些情况下由于数值误差，本应正定的矩阵可能变得半正定或不定。
该函数通过调整特征值来恢复正定性。
"""
function make_positive_definite(A::AbstractMatrix; epsilon::Float64=1e-6)
    A_sym = Symmetric(Matrix{Float64}(A))
    
    # 特征值分解
    eigen_decomp = eigen(A_sym)
    eigenvalues = eigen_decomp.values
    eigenvectors = eigen_decomp.vectors
    
    # 将负特征值替换为epsilon
    eigenvalues_fixed = max.(eigenvalues, epsilon)
    
    # 重构矩阵
    return eigenvectors * Diagonal(eigenvalues_fixed) * eigenvectors'
end

"""
    quadratic_form(x::AbstractVector, A::AbstractMatrix) -> Float64

计算二次型 x' A x。

# 中文说明
高效计算向量的二次型，常用于似然计算。
"""
function quadratic_form(x::AbstractVector, A::AbstractMatrix)
    return dot(x, A * x)
end

"""
    solve_linear_system(A::AbstractMatrix, b::AbstractVector; 
                       method::Symbol=:auto) -> Vector{Float64}

求解线性系统 Ax = b。

# 参数
- `A`: 系数矩阵
- `b`: 右侧向量
- `method`: 求解方法 (:auto, :cholesky, :lu, :qr)

# 中文说明
根据矩阵特性选择最优的求解方法。
对于正定矩阵使用Cholesky分解，否则使用LU分解。
"""
function solve_linear_system(A::AbstractMatrix, b::AbstractVector;
                            method::Symbol=:auto)
    if method == :auto
        # 尝试Cholesky（对于正定对称矩阵最快）
        try
            C = cholesky(Symmetric(A))
            return C \ b
        catch
            # 回退到LU分解
            return lu(A) \ b
        end
    elseif method == :cholesky
        return cholesky(Symmetric(A)) \ b
    elseif method == :lu
        return lu(A) \ b
    elseif method == :qr
        return qr(A) \ b
    else
        throw(ArgumentError("未知的求解方法: $method"))
    end
end

"""
    numerical_gradient(f::Function, x::AbstractVector; epsilon::Float64=1e-8) -> Vector{Float64}

使用有限差分计算函数的数值梯度。

# 中文说明
当解析梯度不可用时的备选方法。
使用中心差分以获得更高精度。
"""
function numerical_gradient(f::Function, x::AbstractVector; epsilon::Float64=1e-8)
    n = length(x)
    grad = zeros(n)
    x_plus = copy(x)
    x_minus = copy(x)
    
    for i in 1:n
        x_plus[i] = x[i] + epsilon
        x_minus[i] = x[i] - epsilon
        
        grad[i] = (f(x_plus) - f(x_minus)) / (2 * epsilon)
        
        x_plus[i] = x[i]
        x_minus[i] = x[i]
    end
    
    return grad
end

"""
    weighted_mean(values::AbstractVector, weights::AbstractVector) -> Float64

计算加权平均。

# 中文说明
加权平均，用于EM算法等场景。
"""
function weighted_mean(values::AbstractVector, weights::AbstractVector)
    @assert length(values) == length(weights) "值和权重长度必须相同"
    return sum(values .* weights) / sum(weights)
end

"""
    weighted_variance(values::AbstractVector, weights::AbstractVector) -> Float64

计算加权方差。

# 中文说明
加权方差估计。
"""
function weighted_variance(values::AbstractVector, weights::AbstractVector)
    mu = weighted_mean(values, weights)
    total_weight = sum(weights)
    return sum(weights .* (values .- mu).^2) / total_weight
end

"""
    robust_inverse(A::AbstractMatrix; regularization::Float64=1e-6) -> Matrix{Float64}

鲁棒的矩阵求逆。

# 中文说明
当矩阵接近奇异时添加正则化项以确保可逆性。
"""
function robust_inverse(A::AbstractMatrix; regularization::Float64=1e-6)
    n = size(A, 1)
    try
        return inv(A)
    catch
        # 添加正则化
        A_reg = A + regularization * I(n)
        return inv(A_reg)
    end
end
