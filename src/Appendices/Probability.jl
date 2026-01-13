# =============================================================================
# Probability.jl - 概率与分布工具
# =============================================================================
# 对应章节: Walsh 2nd Ed, Appendix A2
# =============================================================================

"""
    moment_generating_function(μ::Float64, σ::Float64, t::Float64) -> Float64

正态分布的矩生成函数 M(t) = exp(μt + σ²t²/2)。
"""
function moment_generating_function(μ::Float64, σ::Float64, t::Float64)
    return exp(μ * t + σ^2 * t^2 / 2)
end

"""
    characteristic_function_normal(μ::Float64, σ::Float64, t::Float64) -> ComplexF64

正态分布特征函数 φ(t) = exp(iμt - σ²t²/2)。
"""
function characteristic_function_normal(μ::Float64, σ::Float64, t::Float64)
    return exp(im * μ * t - σ^2 * t^2 / 2)
end

"""
    multivariate_normal_dist(μ::Vector{Float64}, Σ::Matrix{Float64}) -> Distribution

创建多元正态分布。
"""
function multivariate_normal_dist(μ::Vector{Float64}, Σ::Matrix{Float64})
    return MvNormal(μ, Σ)
end

"""
    wishart_density(W::AbstractMatrix, V::AbstractMatrix, n::Int) -> Float64

Wishart分布密度（对数）。
"""
function wishart_density(W::AbstractMatrix, V::AbstractMatrix, n::Int)
    p = size(W, 1)
    
    log_det_W = logdet(W)
    log_det_V = logdet(V)
    tr_VinvW = tr(inv(V) * W)
    
    # 常数项
    log_const = -(n * p / 2) * log(2) - (n / 2) * log_det_V - logmvgamma(p, n/2)
    
    log_dens = log_const + ((n - p - 1) / 2) * log_det_W - tr_VinvW / 2
    
    return log_dens
end

"""
多元gamma函数的对数。
"""
function logmvgamma(p::Int, a::Float64)
    result = p * (p - 1) / 4 * log(π)
    for j in 1:p
        result += loggamma(a + (1 - j) / 2)
    end
    return result
end

"""
    conditional_distribution_mvn(μ::Vector, Σ::Matrix, 
                                obs_idx::Vector{Int}, obs_vals::Vector{Float64}) -> NamedTuple

多元正态的条件分布。
"""
function conditional_distribution_mvn(μ::Vector, Σ::Matrix,
                                      obs_idx::Vector{Int}, obs_vals::Vector{Float64})
    all_idx = 1:length(μ)
    pred_idx = setdiff(all_idx, obs_idx)
    
    μ_1 = μ[pred_idx]
    μ_2 = μ[obs_idx]
    
    Σ_11 = Σ[pred_idx, pred_idx]
    Σ_12 = Σ[pred_idx, obs_idx]
    Σ_22 = Σ[obs_idx, obs_idx]
    
    Σ_22_inv = inv(Σ_22)
    
    # 条件均值和方差
    μ_cond = μ_1 + Σ_12 * Σ_22_inv * (obs_vals - μ_2)
    Σ_cond = Σ_11 - Σ_12 * Σ_22_inv * Σ_12'
    
    return (mean=μ_cond, covariance=Σ_cond)
end
