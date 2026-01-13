# =============================================================================
# OLS.jl - 普通最小二乘回归
# =============================================================================
#
# 实现回归分析与诊断
#
# 对应章节: Walsh 2nd Ed, Chapter 3 (Covariance, Regression, and Correlation)
# =============================================================================

"""
    OLSResult

普通最小二乘回归的结果结构体。

# 字段
- `beta::Vector{Float64}`: 回归系数
- `se::Vector{Float64}`: 系数标准误
- `t_values::Vector{Float64}`: t统计量
- `p_values::Vector{Float64}`: p值
- `R2::Float64`: 决定系数
- `R2_adj::Float64`: 调整后的决定系数
- `residuals::Vector{Float64}`: 残差
- `fitted::Vector{Float64}`: 拟合值
- `sigma2::Float64`: 残差方差
- `n::Int`: 样本量
- `p::Int`: 参数数量
"""
struct OLSResult
    beta::Vector{Float64}
    se::Vector{Float64}
    t_values::Vector{Float64}
    p_values::Vector{Float64}
    R2::Float64
    R2_adj::Float64
    residuals::Vector{Float64}
    fitted::Vector{Float64}
    sigma2::Float64
    n::Int
    p::Int
end

function Base.show(io::IO, r::OLSResult)
    print(io, "OLSResult(R²=$(round(r.R2, digits=4)), R²_adj=$(round(r.R2_adj, digits=4)), n=$(r.n))")
end

"""
    bols_with_intercept(X::AbstractMatrix, y::AbstractVector) -> OLSResult

带截距项的矩阵形式OLS回归。

# 模型方程
```math
\\mathbf{y} = \\mathbf{X}\\boldsymbol{\\beta} + \\boldsymbol{\\epsilon}
```

解：
```math
\\hat{\\boldsymbol{\\beta}} = (\\mathbf{X}^T\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{y}
```

# 参数
- `X::AbstractMatrix`: 设计矩阵（不包含截距列）
- `y::AbstractVector`: 响应变量

# 返回值
- `OLSResult`: 包含回归系数、标准误、R²等

# 中文说明
普通最小二乘（OLS）回归是线性模型分析的基础。
本函数自动添加截距项（全1列）并计算完整的回归诊断。

关键输出：
- β̂: 回归系数估计
- R²: 解释方差的比例
- R²_adj: 考虑参数数量的调整R²

# 公式来源
Walsh 2nd Ed, Eq 3.5-3.8

# 示例
```julia
X = randn(100, 2)  # 两个预测变量
y = 1.0 .+ 2.0 * X[:, 1] .+ 0.5 * X[:, 2] .+ 0.5 * randn(100)

result = bols_with_intercept(X, y)
result.beta  # ≈ [1.0, 2.0, 0.5]
result.R2    # 高R²
```

参见: [`weighted_ols`](@ref)
"""
function bols_with_intercept(X::AbstractMatrix, y::AbstractVector)
    n = length(y)
    
    if size(X, 1) != n
        throw(ArgumentError("X的行数必须与y的长度相同"))
    end
    
    # 添加截距列
    X_aug = hcat(ones(n), X)
    p = size(X_aug, 2)
    
    if n <= p
        throw(ArgumentError("样本量必须大于参数数量"))
    end
    
    # 正规方程求解
    XtX = X_aug' * X_aug
    Xty = X_aug' * y
    
    # 使用Cholesky分解求解（更稳定）
    beta = XtX \ Xty
    
    # 拟合值和残差
    fitted = X_aug * beta
    residuals = y .- fitted
    
    # 残差方差
    df_resid = n - p
    sigma2 = sum(residuals.^2) / df_resid
    
    # 系数协方差矩阵
    cov_beta = sigma2 * inv(XtX)
    se = sqrt.(diag(cov_beta))
    
    # t统计量和p值
    t_values = beta ./ se
    p_values = [2 * ccdf(TDist(df_resid), abs(t)) for t in t_values]
    
    # R² and 调整后R²
    TSS = sum((y .- mean(y)).^2)
    RSS = sum(residuals.^2)
    R2 = 1 - RSS/TSS
    R2_adj = 1 - (RSS/df_resid) / (TSS/(n-1))
    
    return OLSResult(beta, se, t_values, p_values, R2, R2_adj, 
                     residuals, fitted, sigma2, n, p)
end

"""
    simple_regression(x::AbstractVector, y::AbstractVector) -> NamedTuple

简单线性回归 y = a + bx。

# 返回值
- `a`: 截距
- `b`: 斜率
- `r`: 相关系数
- `r2`: 决定系数

# 中文说明
简单线性回归用于分析两个变量间的线性关系。
斜率 b = Cov(x,y) / Var(x)

# 示例
```julia
result = simple_regression(height, weight)
result.b  # 身高每增加1单位，体重增加b单位
```
"""
function simple_regression(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x和y长度必须相同"))
    end
    
    # 过滤有效观测
    valid_idx = Int[]
    for i in 1:length(x)
        if !ismissing(x[i]) && !ismissing(y[i]) &&
           !(x[i] isa Number && isnan(x[i])) && !(y[i] isa Number && isnan(y[i]))
            push!(valid_idx, i)
        end
    end
    
    n = length(valid_idx)
    if n < 3
        throw(ArgumentError("有效样本量不足"))
    end
    
    x_valid = Float64[x[i] for i in valid_idx]
    y_valid = Float64[y[i] for i in valid_idx]
    
    x_mean = mean(x_valid)
    y_mean = mean(y_valid)
    
    # 计算协方差和方差
    cov_xy = sum((x_valid .- x_mean) .* (y_valid .- y_mean)) / (n - 1)
    var_x = sum((x_valid .- x_mean).^2) / (n - 1)
    var_y = sum((y_valid .- y_mean).^2) / (n - 1)
    
    if var_x ≈ 0
        throw(ArgumentError("x无变异，无法进行回归"))
    end
    
    # 斜率和截距
    b = cov_xy / var_x
    a = y_mean - b * x_mean
    
    # 相关系数
    r = cov_xy / sqrt(var_x * var_y)
    r2 = r^2
    
    # 标准误
    y_pred = a .+ b .* x_valid
    residuals = y_valid .- y_pred
    MSE = sum(residuals.^2) / (n - 2)
    se_b = sqrt(MSE / sum((x_valid .- x_mean).^2))
    se_a = sqrt(MSE * (1/n + x_mean^2 / sum((x_valid .- x_mean).^2)))
    
    return (a=a, b=b, se_a=se_a, se_b=se_b, r=r, r2=r2, n=n)
end

"""
    weighted_ols(X::AbstractMatrix, y::AbstractVector, W::AbstractVector) -> OLSResult

加权最小二乘回归。

# 模型
```math
\\hat{\\boldsymbol{\\beta}} = (\\mathbf{X}^T\\mathbf{W}\\mathbf{X})^{-1}\\mathbf{X}^T\\mathbf{W}\\mathbf{y}
```

# 参数
- `X`: 设计矩阵
- `y`: 响应变量
- `W`: 权重向量（对角权重矩阵的对角线）

# 中文说明
加权最小二乘（WLS）用于处理方差异质性。
权重通常设为方差的倒数：w_i = 1/σ²_i

# 公式来源
Walsh 2nd Ed, Chapter 9
"""
function weighted_ols(X::AbstractMatrix, y::AbstractVector, W::AbstractVector)
    n = length(y)
    
    if length(W) != n || size(X, 1) != n
        throw(ArgumentError("维度不匹配"))
    end
    
    # 添加截距
    X_aug = hcat(ones(n), X)
    p = size(X_aug, 2)
    
    # 构建权重矩阵
    W_diag = Diagonal(W)
    
    # 加权正规方程
    XtWX = X_aug' * W_diag * X_aug
    XtWy = X_aug' * W_diag * y
    
    beta = XtWX \ XtWy
    
    # 拟合值和残差
    fitted = X_aug * beta
    residuals = y .- fitted
    
    # 加权残差方差
    df_resid = n - p
    sigma2 = sum(W .* residuals.^2) / df_resid
    
    # 系数协方差
    cov_beta = sigma2 * inv(XtWX)
    se = sqrt.(diag(cov_beta))
    
    # t统计量
    t_values = beta ./ se
    p_values = [2 * ccdf(TDist(df_resid), abs(t)) for t in t_values]
    
    # R²（加权）
    y_mean_w = sum(W .* y) / sum(W)
    TSS_w = sum(W .* (y .- y_mean_w).^2)
    RSS_w = sum(W .* residuals.^2)
    R2 = 1 - RSS_w/TSS_w
    R2_adj = 1 - (RSS_w/df_resid) / (TSS_w/(n-1))
    
    return OLSResult(beta, se, t_values, p_values, R2, R2_adj,
                     residuals, fitted, sigma2, n, p)
end

"""
    parent_offspring_regression(parent::AbstractVector, offspring::AbstractVector) -> NamedTuple

亲子回归估计遗传力。

# 模型
```math
h^2 = 2b_{OP} \\quad \\text{(单亲)}
```
```math
h^2 = b_{MP} \\quad \\text{(双亲均值)}
```

# 返回值
- `b`: 回归斜率
- `h2`: 遗传力估计
- `se_h2`: 遗传力标准误

# 中文说明
亲子回归是估计狭义遗传力的经典方法。

单亲回归：h² = 2 × b（因亲子共享1/2基因）
双亲均值回归：h² = b

# 公式来源
Walsh 2nd Ed, Chapter 7

# 示例
```julia
result = parent_offspring_regression(parent_height, offspring_height)
result.h2  # 遗传力估计
```
"""
function parent_offspring_regression(parent::AbstractVector, offspring::AbstractVector;
                                    single_parent::Bool=true)
    result = simple_regression(parent, offspring)
    
    if single_parent
        # 单亲回归
        h2 = 2 * result.b
        se_h2 = 2 * result.se_b
    else
        # 双亲均值回归
        h2 = result.b
        se_h2 = result.se_b
    end
    
    return (b=result.b, h2=h2, se_h2=se_h2, r=result.r, n=result.n)
end

"""
    partial_correlation_recursion(r12::Float64, r13::Float64, r23::Float64) -> Float64

计算偏相关系数 r₁₂.₃。

# 数学定义
```math
r_{12.3} = \\frac{r_{12} - r_{13}r_{23}}{\\sqrt{(1-r_{13}^2)(1-r_{23}^2)}}
```

# 中文说明
偏相关系数衡量排除第三变量影响后，两个变量间的相关。
r₁₂.₃ 是控制变量3后，变量1和2的相关。

# 公式来源
Walsh 2nd Ed, Eq 3.12

# 示例
```julia
r_xy_z = partial_correlation_recursion(0.7, 0.5, 0.6)  # x和y在控制z后的相关
```
"""
function partial_correlation_recursion(r12::Float64, r13::Float64, r23::Float64)
    denom = sqrt((1 - r13^2) * (1 - r23^2))
    
    if denom ≈ 0
        return NaN
    end
    
    return (r12 - r13 * r23) / denom
end
