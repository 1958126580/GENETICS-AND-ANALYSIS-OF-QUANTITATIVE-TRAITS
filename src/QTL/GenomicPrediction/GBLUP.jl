# =============================================================================
# GBLUP.jl - 基因组BLUP
# =============================================================================
# 对应章节: Walsh 2nd Ed, Chapter 21
# =============================================================================

"""
    gblup_prediction_solver(y::AbstractVector, K::AbstractMatrix;
                           X::AbstractMatrix=ones(length(y),1)) -> PredictionResult

GBLUP基因组预测。

# 模型
y = Xβ + g + ε, g ~ N(0, Gσ²_g)
"""
function gblup_prediction_solver(y::AbstractVector, K::AbstractMatrix;
                                 X::AbstractMatrix=ones(length(y), 1))
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    K_valid = K[valid_idx, valid_idx]
    X_valid = X[valid_idx, :]
    n = length(y_valid)
    
    # 估计方差组分
    vc = reml_ai_algorithm(y_valid, X_valid, K_valid)
    lambda = vc.sigma2_e / vc.sigma2_a
    
    # 求解BLUP
    Z = Matrix{Float64}(I, n, n)
    result = mixed_model_equations(y_valid, X_valid, Z, K_valid, lambda)
    
    gebv = result.u
    
    # 准确度（训练集）
    accuracy = cor(gebv, y_valid .- X_valid * result.beta)
    
    return PredictionResult(gebv, ["ind_$i" for i in valid_idx];
                           accuracy=accuracy, method=:GBLUP)
end

"""
    rr_blup_marker_effects(y::AbstractVector, Z::AbstractMatrix;
                          lambda::Float64=1.0) -> Vector{Float64}

Ridge回归BLUP估计标记效应。

# 模型
y = Zα + ε, E[αα'] = Iσ²_α

α̂ = (Z'Z + λI)⁻¹Z'y
"""
function rr_blup_marker_effects(y::AbstractVector, Z::AbstractMatrix;
                                lambda::Float64=1.0)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    y_valid = Float64[y[i] for i in valid_idx]
    Z_valid = Float64.(Z[valid_idx, :])
    
    n, m = size(Z_valid)
    
    # 中心化
    y_c = y_valid .- mean(y_valid)
    Z_c = Z_valid .- mean(Z_valid, dims=1)
    
    # Ridge回归
    ZtZ = Z_c' * Z_c
    Zty = Z_c' * y_c
    
    alpha = (ZtZ + lambda * I(m)) \ Zty
    
    return alpha
end

"""
    cross_validation_scheme(y::AbstractVector, K::AbstractMatrix;
                           folds::Int=5) -> NamedTuple

K折交叉验证评估预测准确度。
"""
function cross_validation_scheme(y::AbstractVector, K::AbstractMatrix;
                                 folds::Int=5)
    valid_idx = [i for i in 1:length(y) if !ismissing(y[i])]
    n = length(valid_idx)
    
    # 随机分组
    Random.seed!(DEFAULT_RANDOM_SEED)
    fold_assignment = rand(1:folds, n)
    
    predictions = zeros(n)
    observed = Float64[y[valid_idx[i]] for i in 1:n]
    
    for k in 1:folds
        train_idx = findall(f -> f != k, fold_assignment)
        test_idx = findall(f -> f == k, fold_assignment)
        
        y_train = observed[train_idx]
        K_train = K[valid_idx[train_idx], valid_idx[train_idx]]
        K_test_train = K[valid_idx[test_idx], valid_idx[train_idx]]
        
        # 训练
        result = gblup_prediction_solver(y_train, K_train)
        
        # 预测
        gebv_train = result.gebv
        pred_test = K_test_train * (K_train \ gebv_train)
        
        predictions[test_idx] = pred_test
    end
    
    accuracy = cor(predictions, observed)
    
    return (accuracy=accuracy, predictions=predictions, observed=observed)
end
