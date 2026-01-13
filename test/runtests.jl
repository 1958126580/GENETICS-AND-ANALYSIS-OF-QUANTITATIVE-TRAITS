# =============================================================================
# QuantitativeGenetics.jl 完整测试套件
# =============================================================================

using Test
using QuantitativeGenetics
using LinearAlgebra
using Statistics
using Random
using Distributions

Random.seed!(42)

@testset "QuantitativeGenetics.jl Complete Test Suite" begin

    # =========================================================================
    # 类型测试
    # =========================================================================
    @testset "Types" begin
        @testset "GeneticsTypes" begin
            # Allele
            a1 = Allele('A')
            a2 = Allele('a')
            @test a1.id == 'A'
            @test a2.id == 'a'
            
            # Genotype
            g_homo = Genotype(a1, a1)
            g_het = Genotype(a1, a2)
            
            @test is_homozygous(g_homo)
            @test !is_heterozygous(g_homo)
            @test !is_homozygous(g_het)
            @test is_heterozygous(g_het)
            
            # 数值编码
            @test genotype_to_numeric(Genotype(Allele(0), Allele(0))) == 0
            @test genotype_to_numeric(Genotype(Allele(0), Allele(1))) == 1
            @test genotype_to_numeric(Genotype(Allele(1), Allele(1))) == 2
        end
        
        @testset "MapTypes" begin
            m1 = Marker("SNP1", "Chr1", 10.5)
            m2 = Marker("SNP2", "Chr1", 25.0, ["A", "G"])
            
            @test m1.name == "SNP1"
            @test m1.chromosome == "Chr1"
            @test m1.position == 10.5
            @test is_biallelic(m2)
            @test n_alleles(m2) == 2
            
            # GeneticMap
            gmap = GeneticMap("test", [m1, m2])
            @test n_markers(gmap) == 2
            @test n_chromosomes(gmap) == 1
            @test marker_distance(gmap, 1, 2) ≈ 14.5
        end
        
        @testset "PopulationTypes" begin
            # TraitData
            td = TraitData([1.0, 2.0, missing, 4.0])
            @test n_observed(td) == 3
            @test trait_mean(td) ≈ 7/3 atol=0.01
            
            # Population
            geno = rand(0:2, 10, 5)
            traits = Dict("height" => Vector{Union{Float64,Missing}}(randn(10)))
            pop = Population("test", ["ind_$i" for i in 1:10], 
                            ["M$i" for i in 1:5], Float64.(geno), traits, :F2)
            
            @test n_individuals(pop) == 10
            @test n_markers(pop) == 5
            @test n_traits(pop) == 1
        end
        
        @testset "ModelTypes" begin
            add_model = AdditiveModel(0.5, 1.0)
            @test add_model.mu == 0.5
            @test add_model.a == 1.0
            
            dom_model = DominanceModel(0.0, 1.0, 0.5)
            @test dominance_ratio(dom_model) == 0.5
        end
        
        @testset "RelationshipTypes" begin
            ped = Pedigree(
                ["A", "B", "C", "D"],
                [nothing, nothing, "A", "A"],
                [nothing, nothing, "B", "B"]
            )
            @test n_individuals(ped) == 4
            @test n_founders(ped) == 2
            @test is_founder(ped, "A")
            @test !is_founder(ped, "C")
        end
        
        @testset "ResultTypes" begin
            hwe_res = HWEResult([25,50,25], [25.0,50.0,25.0], 0.0, 1.0, :chi_square, (0.5,0.5))
            @test is_in_hwe(hwe_res)
            
            vc = VarianceComponents(0.5, 0.5)
            @test vc.heritability ≈ 0.5 atol=0.01
        end
    end

    # =========================================================================
    # 频率与HWE测试
    # =========================================================================
    @testset "Allele Frequencies and HWE" begin
        @testset "Frequency Calculation" begin
            result = allele_genotype_freqs([25, 50, 25])
            @test result.p ≈ 0.5 atol=0.001
            @test result.q ≈ 0.5 atol=0.001
            @test sum(result.geno_freqs) ≈ 1.0
            
            # 从基因型矩阵
            genos = [0, 1, 1, 2, 0, 1, 0, 2]
            result2 = allele_freqs_from_matrix(Float64.(genos))
            @test 0 < result2.maf < 0.5
        end
        
        @testset "HWE Expected" begin
            hwe = hwe_expected_freqs(0.3, 100)
            @test hwe.exp_freqs[1] ≈ 0.09 atol=0.001  # p²
            @test hwe.exp_freqs[2] ≈ 0.42 atol=0.001  # 2pq
            @test hwe.exp_freqs[3] ≈ 0.49 atol=0.001  # q²
        end
        
        @testset "HWE Tests" begin
            # 符合HWE的数据
            result_hwe = hwe_chi_square_test([25, 50, 25])
            @test result_hwe.pvalue > 0.05
            
            # 偏离HWE的数据
            result_dev = hwe_chi_square_test([40, 20, 40])
            @test result_dev.pvalue < 0.05
            
            # 精确检验
            exact_result = hwe_exact_test([10, 20, 10])
            @test exact_result.pvalue > 0.05
        end
        
        @testset "Heterozygosity" begin
            h = heterozygosity(0.5)
            @test h.H_exp ≈ 0.5 atol=0.001
            @test h.H_max == 0.5
            
            h_multi = heterozygosity_multiallelic([0.5, 0.3, 0.2])
            @test h_multi ≈ 1 - (0.25 + 0.09 + 0.04) atol=0.001
        end
    end

    # =========================================================================
    # Fisher分解与育种值
    # =========================================================================
    @testset "Fisher Decomposition and Breeding Values" begin
        @testset "Fisher Partition" begin
            result = fisher_decomposition(0.5, 0.5, 1.0, 0.5)
            @test result.alpha ≈ 1.0 atol=0.01
            @test result.V_A > 0
            @test result.V_D >= 0
            @test result.V_G ≈ result.V_A + result.V_D atol=0.001
        end
        
        @testset "Breeding Values" begin
            bv = breeding_value_partition(0.5, 0.5, 1.0)
            @test bv.A_AA ≈ 1.0 atol=0.01  # 2qα
            @test bv.A_Aa ≈ 0.0 atol=0.01  # (q-p)α
            @test bv.A_aa ≈ -1.0 atol=0.01 # -2pα
        end
        
        @testset "Dominance Deviation" begin
            dd = dominance_deviation_calc(0.5, 0.5, 0.5)
            @test dd.D_Aa > 0  # 杂合子正偏差
        end
    end

    # =========================================================================
    # 回归分析
    # =========================================================================
    @testset "Regression Analysis" begin
        @testset "Simple Regression" begin
            x = collect(1.0:10.0)
            y = 2.0 .+ 3.0 .* x .+ 0.1 .* randn(10)
            
            result = simple_regression(x, y)
            @test result.b ≈ 3.0 atol=0.2
            @test result.a ≈ 2.0 atol=0.5
            @test result.r2 > 0.99
        end
        
        @testset "OLS with Intercept" begin
            n = 100
            X = randn(n, 2)
            y = 1.0 .+ 2.0 .* X[:, 1] .+ 0.5 .* X[:, 2] .+ 0.3 .* randn(n)
            
            result = bols_with_intercept(X, y)
            @test result.beta[1] ≈ 1.0 atol=0.3
            @test result.beta[2] ≈ 2.0 atol=0.3
            @test result.beta[3] ≈ 0.5 atol=0.3
            @test result.R2 > 0.9
        end
        
        @testset "Parent-Offspring Regression" begin
            parent = randn(50) .+ 10.0
            offspring = 0.6 .* parent .+ 0.4 .* randn(50) .* std(parent) .+ 4.0
            
            result = parent_offspring_regression(parent, offspring)
            @test 0.3 < result.h2 < 0.9  # h² = 2b for single parent
        end
    end

    # =========================================================================
    # 连锁不平衡
    # =========================================================================
    @testset "Linkage Disequilibrium" begin
        @testset "LD Metrics" begin
            result = linkage_disequilibrium_metrics(0.3, 0.5, 0.5)
            @test result.D ≈ 0.05 atol=0.01
            @test -1 <= result.D_prime <= 1
            @test 0 <= result.r_squared <= 1
        end
        
        @testset "LD from Genotypes" begin
            n = 100
            g1 = rand(0:2, n)
            g2 = copy(g1) .+ rand([-1, 0, 1], n)
            g2 = clamp.(g2, 0, 2)
            
            result = ld_from_genotypes(Float64.(g1), Float64.(g2))
            @test result.r_squared > 0.2  # 应该有相关性
        end
    end

    # =========================================================================
    # 映射函数
    # =========================================================================
    @testset "Map Functions" begin
        @testset "Haldane" begin
            d = haldane_mapping(0.1)
            r = inverse_haldane(d)
            @test r ≈ 0.1 atol=0.001
            
            @test haldane_mapping(0.0) ≈ 0.0 atol=0.01
        end
        
        @testset "Kosambi" begin
            d_k = kosambi_mapping(0.1)
            r_k = inverse_kosambi(d_k)
            @test r_k ≈ 0.1 atol=0.001
        end
        
        @testset "Comparison" begin
            r = 0.2
            d_h = haldane_mapping(r)
            d_k = kosambi_mapping(r)
            @test d_k < d_h  # Kosambi假设正干扰，距离更短
        end
    end

    # =========================================================================
    # 混合模型
    # =========================================================================
    @testset "Mixed Models" begin
        @testset "MME" begin
            n = 30
            y = randn(n)
            X = ones(n, 1)
            Z = Matrix{Float64}(I, n, n)
            K = 0.1 * ones(n, n) + 0.9 * I(n)
            
            result = mixed_model_equations(y, X, Z, K, 1.0)
            @test length(result.beta) == 1
            @test length(result.u) == n
        end
        
        @testset "REML" begin
            n = 50
            y = randn(n)
            X = ones(n, 1)
            K = Matrix{Float64}(I, n, n) + 0.1 * ones(n, n)
            
            vc = reml_ai_algorithm(y, X, K; max_iter=30)
            @test vc.sigma2_a >= 0
            @test vc.sigma2_e >= 0
            @test 0 <= vc.heritability <= 1
        end
    end

    # =========================================================================
    # QTL作图
    # =========================================================================
    @testset "QTL Mapping" begin
        @testset "Haley-Knott Scan" begin
            n, m = 100, 10
            geno = rand(0:2, n, m)
            qtl_idx = 5
            y = Float64.(geno[:, qtl_idx]) .+ 0.5 .* randn(n)
            
            traits = Dict("trait" => Vector{Union{Float64,Missing}}(y))
            pop = Population("test", ["ind_$i" for i in 1:n],
                            ["M$i" for i in 1:m], Float64.(geno), traits, :F2)
            
            result = haley_knott_scan(pop, :trait)
            @test length(result.lod) == m
            @test argmax(result.lod) == qtl_idx
        end
    end

    # =========================================================================
    # GWAS
    # =========================================================================
    @testset "GWAS" begin
        @testset "Single Marker Scan" begin
            n, m = 100, 20
            y = randn(n)
            geno = rand(0:2, n, m)
            
            result = gwas_single_marker_scan(y, Float64.(geno))
            @test length(result.pvalues) == m
            @test all(0 .<= result.pvalues .<= 1)
            @test result.lambda_gc > 0
        end
        
        @testset "Genomic Control" begin
            pvals = rand(100)
            lambda = genomic_control_lambda(pvals)
            @test 0.5 < lambda < 2.0  # 应该接近1
        end
        
        @testset "FDR Correction" begin
            pvals = [0.001, 0.01, 0.05, 0.1, 0.5]
            result = fdr_correction(pvals; alpha=0.1)
            @test length(result.q_values) == 5
        end
    end

    # =========================================================================
    # 基因组预测
    # =========================================================================
    @testset "Genomic Prediction" begin
        @testset "GBLUP" begin
            n = 50
            y = randn(n)
            K = Matrix{Float64}(I, n, n) + 0.2 * randn(n, n)
            K = (K + K') / 2 + 0.5 * I(n)
            
            result = gblup_prediction_solver(y, K)
            @test length(result.gebv) == n
            @test !isnan(result.accuracy)
        end
        
        @testset "RR-BLUP" begin
            n, m = 50, 100
            y = randn(n)
            Z = randn(n, m)
            
            effects = rr_blup_marker_effects(y, Z; lambda=10.0)
            @test length(effects) == m
        end
    end

    # =========================================================================
    # 品系杂交
    # =========================================================================
    @testset "Line Crosses" begin
        @testset "Generation Mean Analysis" begin
            P1, P2 = 100.0, 80.0
            F1 = 95.0  # 显性偏向P1
            F2 = 92.0
            BC1 = 97.0
            BC2 = 87.0
            
            result = generation_mean_analysis(P1, P2, F1, F2, BC1, BC2)
            @test result.m ≈ 90.0 atol=5.0  # 中亲值
            @test result.a > 0  # 加性效应
        end
        
        @testset "Castle-Wright" begin
            n_genes = castle_wright_n_factors(100.0, 80.0, 50.0, 10.0)
            @test n_genes > 0
        end
    end

    # =========================================================================
    # 近交与杂种优势
    # =========================================================================
    @testset "Inbreeding and Heterosis" begin
        @testset "Inbreeding Depression" begin
            mu = inbreeding_depression_linear(100.0, 0.25, 20.0)
            @test mu ≈ 95.0 atol=0.01
        end
        
        @testset "Heterosis" begin
            result = dominance_vs_overdominance_heterosis(100.0, 80.0, 95.0)
            @test result.mp_heterosis ≈ 5.0 atol=0.01
            @test result.dom_ratio < 1.0  # 部分显性
        end
    end

    # =========================================================================
    # 数值工具
    # =========================================================================
    @testset "Numerical Utils" begin
        @testset "Safe Functions" begin
            @test safe_log(0.0) ≈ log(MIN_PROBABILITY)
            @test safe_log(1.0) ≈ 0.0
            
            @test safe_exp(0.0) ≈ 1.0
            @test safe_exp(-1000.0) ≈ 0.0 atol=1e-10
        end
        
        @testset "Logit/Expit" begin
            @test logit(0.5) ≈ 0.0 atol=1e-10
            @test expit(0.0) ≈ 0.5 atol=1e-10
            
            for p in [0.1, 0.3, 0.5, 0.7, 0.9]
                @test expit(logit(p)) ≈ p atol=1e-10
            end
        end
        
        @testset "LogSumExp" begin
            x = [-1000.0, -1001.0, -1002.0]
            result = logsumexp(x)
            @test result > -1000.0
            @test !isinf(result)
        end
    end

    # =========================================================================
    # 矩阵代数
    # =========================================================================
    @testset "Matrix Algebra" begin
        @testset "Kronecker Product" begin
            A = [1.0 2.0; 3.0 4.0]
            B = [1.0 0.0; 0.0 1.0]
            K = kronecker_product(A, B)
            @test size(K) == (4, 4)
        end
        
        @testset "Generalized Inverse" begin
            A = randn(3, 3)
            A_pinv = generalized_inverse(A)
            @test size(A_pinv) == (3, 3)
        end
    end

    # =========================================================================
    # 亲缘关系矩阵
    # =========================================================================
    @testset "Relationship Matrices" begin
        @testset "Pedigree-based A matrix" begin
            ped = Pedigree(
                ["S", "D", "O1", "O2"],
                [nothing, nothing, "S", "S"],
                [nothing, nothing, "D", "D"]
            )
            
            A = additive_relationship_matrix(ped)
            @test A.matrix[1, 1] ≈ 1.0 atol=0.01  # 自身
            @test A.matrix[3, 4] ≈ 0.5 atol=0.01  # 全同胞
        end
        
        @testset "GRM" begin
            n, m = 20, 50
            geno = rand(0:2, n, m)
            p = [mean(geno[:, j]) / 2 for j in 1:m]
            
            G = grm_vanraden_yang(Float64.(geno), p)
            @test size(G.matrix) == (n, n)
            @test issymmetric(G.matrix)
        end
    end

end  # End main testset

println("\n✅ All tests completed!")
