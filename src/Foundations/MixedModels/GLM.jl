# =============================================================================
# GLM.jl - 广义线性模型
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 9
# =============================================================================

"""
    general_linear_model(y::AbstractVector, X::AbstractMatrix;
                        W::Union{AbstractVector, Nothing}=nothing) -> OLSResult

广义线性模型（加权最小二乘）。

# 模型
y = Xβ + ε, ε ~ N(0, W⁻¹σ²)

# 中文说明
加权最小二乘用于处理异方差数据。
"""
function general_linear_model(y::AbstractVector, X::AbstractMatrix;
                              W::Union{AbstractVector, Nothing}=nothing)
    if W === nothing
        return bols_with_intercept(X, y)
    else
        return weighted_ols(X, y, W)
    end
end
