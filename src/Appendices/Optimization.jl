# =============================================================================
# Optimization.jl - 优化算法
# =============================================================================
# 对应章节: Walsh 2nd Ed, Appendix A3
# =============================================================================

"""
    newton_raphson_update(f::Function, grad_f::Function, hess_f::Function,
                         x0::Vector{Float64}; max_iter::Int=100, tol::Float64=1e-6) -> NamedTuple

Newton-Raphson优化。

x_{k+1} = x_k - H(x_k)⁻¹ ∇f(x_k)
"""
function newton_raphson_update(f::Function, grad_f::Function, hess_f::Function,
                               x0::Vector{Float64}; max_iter::Int=100, tol::Float64=1e-6)
    x = copy(x0)
    converged = false
    
    for i in 1:max_iter
        g = grad_f(x)
        H = hess_f(x)
        
        # 更新
        delta = -H \ g
        x_new = x + delta
        
        if norm(delta) < tol
            converged = true
            x = x_new
            break
        end
        
        x = x_new
    end
    
    return (solution=x, value=f(x), converged=converged)
end

"""
    em_algorithm_template(E_step::Function, M_step::Function, params0;
                         max_iter::Int=100, tol::Float64=1e-6) -> NamedTuple

EM算法模板。
"""
function em_algorithm_template(E_step::Function, M_step::Function, params0;
                               max_iter::Int=100, tol::Float64=1e-6)
    params = params0
    log_lik_prev = -Inf
    converged = false
    
    for i in 1:max_iter
        # E步
        expectations = E_step(params)
        
        # M步
        params_new, log_lik = M_step(expectations)
        
        # 检查收敛
        if abs(log_lik - log_lik_prev) < tol
            converged = true
            params = params_new
            break
        end
        
        params = params_new
        log_lik_prev = log_lik
    end
    
    return (params=params, converged=converged)
end

"""
    fisher_scoring_iteration(score::Function, fisher_info::Function,
                            theta0::Vector{Float64}; max_iter::Int=50) -> Vector{Float64}

Fisher Scoring算法。

θ_{k+1} = θ_k + I(θ_k)⁻¹ × U(θ_k)
"""
function fisher_scoring_iteration(score::Function, fisher_info::Function,
                                  theta0::Vector{Float64}; max_iter::Int=50, tol::Float64=1e-6)
    theta = copy(theta0)
    
    for i in 1:max_iter
        U = score(theta)
        I_mat = fisher_info(theta)
        
        delta = I_mat \ U
        theta_new = theta + delta
        
        if norm(delta) < tol
            return theta_new
        end
        
        theta = theta_new
    end
    
    return theta
end

"""
    gradient_descent_basic(f::Function, grad::Function, x0::Vector{Float64};
                          lr::Float64=0.01, max_iter::Int=1000) -> Vector{Float64}

基础梯度下降。
"""
function gradient_descent_basic(f::Function, grad::Function, x0::Vector{Float64};
                                lr::Float64=0.01, max_iter::Int=1000, tol::Float64=1e-6)
    x = copy(x0)
    
    for i in 1:max_iter
        g = grad(x)
        x_new = x - lr * g
        
        if norm(g) < tol
            return x_new
        end
        
        x = x_new
    end
    
    return x
end
