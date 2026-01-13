# =============================================================================
# MME.jl - Henderson混合模型方程
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 10
# =============================================================================

"""
    mixed_model_equations(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix,
                         K::AbstractMatrix, lambda::Float64) -> NamedTuple

Henderson混合模型方程求解。

# 模型
y = Xβ + Zu + ε

# MME系统
[X'X    X'Z  ] [β̂]   [X'y]
[Z'X  Z'Z+λK⁻¹] [û] = [Z'y]

其中 λ = σ²_e / σ²_u

# 返回
- beta: 固定效应BLUE
- u: 随机效应BLUP
- V: 方差协方差矩阵
"""
function mixed_model_equations(y::AbstractVector, X::AbstractMatrix, 
                               Z::AbstractMatrix, K::AbstractMatrix,
                               lambda::Float64)
    n = length(y)
    p = size(X, 2)
    q = size(Z, 2)
    
    @assert size(X, 1) == n "X行数与y不匹配"
    @assert size(Z, 1) == n "Z行数与y不匹配"
    @assert size(K) == (q, q) "K矩阵维度与Z不匹配"
    
    # K的逆（加正则化）
    K_inv = inv(K + 1e-6 * I(q))
    
    # 构建MME系数矩阵
    C11 = X' * X
    C12 = X' * Z
    C21 = Z' * X
    C22 = Z' * Z + lambda * K_inv
    
    LHS = [C11 C12; C21 C22]
    RHS = [X' * y; Z' * y]
    
    # 求解
    solutions = LHS \ RHS
    
    beta = solutions[1:p]
    u = solutions[(p+1):end]
    
    return (beta=beta, u=u, LHS=LHS, RHS=RHS)
end

"""
    blue_blup_decomposition(mme_solution::NamedTuple) -> NamedTuple

提取MME解中的BLUE和BLUP。
"""
function blue_blup_decomposition(mme_solution::NamedTuple)
    return (BLUE=mme_solution.beta, BLUP=mme_solution.u)
end

"""
    animal_model_scan(y::AbstractVector, A::AbstractMatrix, 
                     X::AbstractMatrix=ones(length(y),1);
                     sigma2_a::Float64=1.0, sigma2_e::Float64=1.0) -> NamedTuple

动物模型BLUP评估。

# 模型
y = Xb + Za + e
a ~ N(0, Aσ²_a)

# 中文说明
动物模型是育种值预测的标准方法。
"""
function animal_model_scan(y::AbstractVector, A::AbstractMatrix;
                          X::AbstractMatrix=ones(length(y), 1),
                          sigma2_a::Float64=1.0, sigma2_e::Float64=1.0)
    n = length(y)
    lambda = sigma2_e / sigma2_a
    
    # Z = I (每个个体有自己的育种值)
    Z = Matrix{Float64}(I, n, n)
    
    result = mixed_model_equations(y, X, Z, A, lambda)
    
    return (fixed_effects=result.beta, 
            breeding_values=result.u,
            reliability=1.0 .- diag(inv(result.LHS)[(length(result.beta)+1):end, 
                                                     (length(result.beta)+1):end]) .* sigma2_a)
end
