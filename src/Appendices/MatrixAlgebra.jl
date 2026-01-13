# =============================================================================
# MatrixAlgebra.jl - 矩阵代数工具
# =============================================================================
# 对应章节: Walsh 2nd Ed, Appendix A1
# =============================================================================

"""
    kronecker_product(A::AbstractMatrix, B::AbstractMatrix) -> Matrix

Kronecker积 A ⊗ B。
"""
function kronecker_product(A::AbstractMatrix, B::AbstractMatrix)
    m1, n1 = size(A)
    m2, n2 = size(B)
    
    result = zeros(m1 * m2, n1 * n2)
    for i in 1:m1
        for j in 1:n1
            result[(i-1)*m2+1:i*m2, (j-1)*n2+1:j*n2] = A[i, j] * B
        end
    end
    
    return result
end

"""
    hadamard_product(A::AbstractMatrix, B::AbstractMatrix) -> Matrix

Hadamard积（逐元素乘积）A ∘ B。
"""
hadamard_product(A::AbstractMatrix, B::AbstractMatrix) = A .* B

"""
    vec_operator(A::AbstractMatrix) -> Vector

矩阵向量化算子 vec(A)。
"""
vec_operator(A::AbstractMatrix) = vec(A)

"""
    vech_operator(A::AbstractMatrix) -> Vector

下三角向量化算子 vech(A)。
"""
function vech_operator(A::AbstractMatrix)
    n = size(A, 1)
    result = Float64[]
    for j in 1:n
        for i in j:n
            push!(result, A[i, j])
        end
    end
    return result
end

"""
    duplication_matrix(n::Int) -> Matrix

复制矩阵 D_n，满足 vec(A) = D_n × vech(A)。
"""
function duplication_matrix(n::Int)
    vech_len = n * (n + 1) ÷ 2
    vec_len = n^2
    D = zeros(vec_len, vech_len)
    
    k = 0
    for j in 1:n
        for i in j:n
            k += 1
            # vec位置
            D[(j-1)*n + i, k] = 1.0
            if i != j
                D[(i-1)*n + j, k] = 1.0
            end
        end
    end
    
    return D
end

"""
    commutation_matrix(m::Int, n::Int) -> Matrix

交换矩阵 K_{m,n}，满足 K × vec(A) = vec(A')。
"""
function commutation_matrix(m::Int, n::Int)
    K = zeros(m * n, m * n)
    for i in 1:m
        for j in 1:n
            # vec(A)中位置 (i,j) 映射到 vec(A')中位置 (j,i)
            from_idx = (j-1)*m + i
            to_idx = (i-1)*n + j
            K[to_idx, from_idx] = 1.0
        end
    end
    return K
end

"""
    generalized_inverse(A::AbstractMatrix; tol::Float64=1e-10) -> Matrix

Moore-Penrose广义逆 A⁺。
"""
function generalized_inverse(A::AbstractMatrix; tol::Float64=1e-10)
    svd_result = svd(A)
    
    # 对奇异值求倒数（小于tol的设为0）
    d_inv = [s > tol ? 1/s : 0.0 for s in svd_result.S]
    
    return svd_result.Vt' * Diagonal(d_inv) * svd_result.U'
end

"""
    matrix_differential_rules(A::AbstractMatrix) -> NamedTuple

返回常用矩阵微分公式参考。
"""
function matrix_differential_rules(A::AbstractMatrix)
    n = size(A, 1)
    return (
        trace_derivative = "d(tr(AX))/dX = A'",
        det_derivative = "d(|A|)/dA = |A|(A⁻¹)'",
        inverse_derivative = "d(A⁻¹)/dA = -A⁻¹ ⊗ (A⁻¹)'"
    )
end
