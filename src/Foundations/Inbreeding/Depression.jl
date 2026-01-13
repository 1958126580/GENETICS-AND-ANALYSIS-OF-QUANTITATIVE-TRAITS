# =============================================================================
# Depression.jl - 近交衰退与杂种优势
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapters 12-13
# =============================================================================

"""
    inbreeding_depression_linear(mu0::Float64, F::Float64, delta::Float64) -> Float64

线性近交衰退模型。

# 公式
μ_F = μ_0 - δF

# 参数
- mu0: 非近交群体均值
- F: 近交系数
- delta: 近交衰退系数
"""
inbreeding_depression_linear(mu0::Float64, F::Float64, delta::Float64) = mu0 - delta * F

"""
    estimate_lethal_equivalents(survival::Vector{Float64}, 
                               F::Vector{Float64}) -> NamedTuple

估计致死当量数B。

# 模型
ln(survival) = A - BF

B是致死当量数（每单位近交导致的适应度下降）。
"""
function estimate_lethal_equivalents(survival::Vector{Float64}, F::Vector{Float64})
    @assert length(survival) == length(F)
    @assert all(s -> s > 0, survival)
    
    log_surv = log.(survival)
    result = simple_regression(F, log_surv)
    
    B = -result.b
    A = exp(result.a)
    
    return (B=B, A=A, r_squared=result.r2)
end

"""
    dominance_vs_overdominance_heterosis(P1::Float64, P2::Float64, 
                                         F1::Float64) -> NamedTuple

判断杂种优势的遗传基础。

# 指标
- mid_parent_heterosis: F1与双亲均值的偏差
- best_parent_heterosis: F1与优势亲本的偏差

# 判断
- d/|a| < 1: 部分显性
- d/|a| ≈ 1: 完全显性
- d/|a| > 1: 超显性
"""
function dominance_vs_overdominance_heterosis(P1::Float64, P2::Float64, F1::Float64)
    mid_parent = (P1 + P2) / 2
    best_parent = max(P1, P2)
    
    mp_heterosis = F1 - mid_parent
    bp_heterosis = F1 - best_parent
    
    # 估计a和d
    a = (P1 - P2) / 2
    d = F1 - mid_parent
    
    dom_ratio = abs(a) > 0 ? abs(d) / abs(a) : Inf
    
    dominance_type = if dom_ratio < 0.5
        :additive
    elseif dom_ratio < 1.0
        :partial_dominance
    elseif dom_ratio ≈ 1.0
        :complete_dominance
    else
        :overdominance
    end
    
    return (mp_heterosis=mp_heterosis, bp_heterosis=bp_heterosis,
            a=a, d=d, dom_ratio=dom_ratio, dominance_type=dominance_type)
end

"""
    outbreeding_depression_model(F1_hyb::Float64, F2_hyb::Float64) -> Float64

估计远交衰退（杂交衰退）。

F2相对于F1的适应度下降。
"""
function outbreeding_depression_model(F1_hyb::Float64, F2_hyb::Float64)
    return (F1_hyb - F2_hyb) / F1_hyb
end
