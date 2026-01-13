# =============================================================================
# MutationalVariance.jl - 突变方差
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 15
# =============================================================================

"""
    mutational_variance_vm_estimation(V_A_t::Vector{Float64}, 
                                     generations::Vector{Int}) -> NamedTuple

从突变累积实验估计突变方差V_m。

# 模型
V_A(t) = V_A(0) + t × V_m
"""
function mutational_variance_vm_estimation(V_A_t::Vector{Float64},
                                           generations::Vector{Int})
    result = simple_regression(Float64.(generations), V_A_t)
    
    V_m = result.b  # 每代增加的方差
    V_A_0 = result.a  # 初始方差
    
    return (V_m=V_m, V_A_0=V_A_0, r_squared=result.r2)
end
