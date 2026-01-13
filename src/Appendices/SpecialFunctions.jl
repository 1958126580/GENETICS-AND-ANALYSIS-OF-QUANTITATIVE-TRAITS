# =============================================================================
# SpecialFunctions.jl - 特殊函数
# =============================================================================
# 对应章节: Walsh 2nd Ed, Appendix A4
# =============================================================================

"""
    beta_function(a::Float64, b::Float64) -> Float64

Beta函数 B(a,b) = Γ(a)Γ(b)/Γ(a+b)。
"""
beta_function(a::Float64, b::Float64) = exp(loggamma(a) + loggamma(b) - loggamma(a + b))

"""
    incomplete_beta(x::Float64, a::Float64, b::Float64) -> Float64

不完全Beta函数 I_x(a,b)。
"""
incomplete_beta(x::Float64, a::Float64, b::Float64) = beta_inc(a, b, x)[1]

"""
    incomplete_gamma(a::Float64, x::Float64) -> Float64

不完全Gamma函数 γ(a,x)/Γ(a)。
"""
incomplete_gamma(a::Float64, x::Float64) = gamma_inc(a, x)[1]

"""
    hypergeometric_1f1(a::Float64, b::Float64, z::Float64) -> Float64

合流超几何函数 ₁F₁(a;b;z) 的级数近似。
"""
function hypergeometric_1f1(a::Float64, b::Float64, z::Float64; max_terms::Int=100)
    result = 1.0
    term = 1.0
    
    for n in 1:max_terms
        term *= (a + n - 1) * z / ((b + n - 1) * n)
        result += term
        
        if abs(term) < 1e-12
            break
        end
    end
    
    return result
end

"""
    hermite_polynomial(n::Int, x::Float64) -> Float64

物理学家的Hermite多项式 Hₙ(x)。
"""
function hermite_polynomial(n::Int, x::Float64)
    if n == 0
        return 1.0
    elseif n == 1
        return 2x
    else
        H_prev2 = 1.0
        H_prev1 = 2x
        H = 0.0
        for k in 2:n
            H = 2x * H_prev1 - 2(k-1) * H_prev2
            H_prev2 = H_prev1
            H_prev1 = H
        end
        return H
    end
end

"""
    laguerre_polynomial(n::Int, x::Float64) -> Float64

Laguerre多项式 Lₙ(x)。
"""
function laguerre_polynomial(n::Int, x::Float64)
    if n == 0
        return 1.0
    elseif n == 1
        return 1 - x
    else
        L_prev2 = 1.0
        L_prev1 = 1 - x
        L = 0.0
        for k in 2:n
            L = ((2k - 1 - x) * L_prev1 - (k - 1) * L_prev2) / k
            L_prev2 = L_prev1
            L_prev1 = L
        end
        return L
    end
end

"""
    bessel_i0(x::Float64) -> Float64

修正Bessel函数 I₀(x) 的近似。
"""
function bessel_i0(x::Float64)
    ax = abs(x)
    if ax < 3.75
        y = (x / 3.75)^2
        return 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
               y*(0.2659732 + y*(0.0360768 + y*0.0045813)))))
    else
        y = 3.75 / ax
        return (exp(ax) / sqrt(ax)) * (0.39894228 + y*(0.01328592 +
               y*(0.00225319 + y*(-0.00157565 + y*(0.00916281 +
               y*(-0.02057706 + y*(0.02635537 + y*(-0.01647633 + y*0.00392377))))))))
    end
end
