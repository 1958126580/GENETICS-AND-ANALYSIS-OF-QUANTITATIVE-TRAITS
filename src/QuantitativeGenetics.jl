# =============================================================================
# QuantitativeGenetics.jl
# =============================================================================
#
# 数量性状遗传学完整实现
# Julia实现 Walsh, Visscher & Lynch (2025) 
# "Genetics and Analysis of Quantitative Traits: Foundations, 2nd Edition"
#
# =============================================================================

module QuantitativeGenetics

# ========== 标准库 ==========
using LinearAlgebra
using Statistics
using Random
using SparseArrays

# ========== 外部依赖 ==========
using Distributions
using SpecialFunctions: loggamma, gamma_inc, beta_inc
using CSV
using DataFrames

# ========== 常量 ==========
include("Constants.jl")

# ========== 类型定义 ==========
include("Types/GeneticsTypes.jl")
include("Types/MapTypes.jl")
include("Types/PopulationTypes.jl")
include("Types/ModelTypes.jl")
include("Types/RelationshipTypes.jl")
include("Types/ResultTypes.jl")

# ========== 工具函数 ==========
include("Utils/NumericalUtils.jl")
include("Utils/Validation.jl")
include("Utils/IO.jl")

# ========== 基础模块 (Chapters 1-14) ==========
# Chapter 1-2: 表型与分布
include("Foundations/Distributions/Phenotypes.jl")
include("Foundations/Distributions/Moments.jl")

# Chapter 3: 回归与通径分析
include("Foundations/Regression/OLS.jl")
include("Foundations/Regression/PathAnalysis.jl")

# Chapter 4: 单位点遗传学
include("Foundations/SingleLoci/Frequencies.jl")
include("Foundations/SingleLoci/HWE.jl")
include("Foundations/SingleLoci/BreedingValue.jl")

# Chapter 5: 上位性
include("Foundations/Epistasis/TwoLocus.jl")

# Chapter 6: 连锁不平衡与环境互作
include("Foundations/LinkageDisequilibrium/LD.jl")
include("Foundations/Environment/GxE.jl")

# Chapter 7-8: 亲缘关系
include("Foundations/Resemblance/IBD.jl")
include("Foundations/Resemblance/GRM.jl")

# Chapter 9-10: 混合模型
include("Foundations/MixedModels/GLM.jl")
include("Foundations/MixedModels/MME.jl")
include("Foundations/MixedModels/REML.jl")

# Chapter 11-14: 品系杂交与尺度
include("Foundations/LineCrosses/GenerationMeans.jl")
include("Foundations/Inbreeding/Depression.jl")
include("Foundations/Scale/Transformations.jl")

# ========== 连锁分析 (Chapter 17) ==========
include("Linkage/MapFunctions.jl")
include("Linkage/MapConstruction.jl")

# ========== QTL作图 (Chapters 15-19) ==========
include("QTL/MajorGenes/Segregation.jl")
include("QTL/Mutation/MutationalVariance.jl")
include("QTL/IntervalMapping/IM.jl")
include("QTL/IntervalMapping/CIM.jl")
include("QTL/Outbred/VCLinkage.jl")

# ========== GWAS与基因组预测 (Chapters 20-21) ==========
include("QTL/GWAS/Association.jl")
include("QTL/GWAS/Visualization.jl")
include("QTL/GenomicPrediction/GBLUP.jl")
include("QTL/GenomicPrediction/BayesAlphabet.jl")

# ========== 数学工具 (Appendices) ==========
include("Appendices/MatrixAlgebra.jl")
include("Appendices/Probability.jl")
include("Appendices/Optimization.jl")
include("Appendices/SpecialFunctions.jl")

# ========== 导出符号 ==========

# 类型
export Allele, Genotype, Haplotype, DiploidGenome
export Marker, GeneticMap, IntervalPosition
export TraitData, Individual, Population, CrossDesign
export AbstractModel, QTLModel, AdditiveModel, DominanceModel, EpistaticModel
export CompositeModel, MixedModel, VarianceComponentsModel
export Pedigree, RelationshipMatrix, IBDProbabilities, Δ_coefficients
export HWEResult, BreedingValueResult, VarianceComponents
export QTLScanResult, GWASResult, PredictionResult, LineCrossResult
export ValidationResult, OLSResult, TraitType

# 常量
export DEFAULT_TOLERANCE, STRICT_TOLERANCE, MAX_ITERATIONS
export MIN_PROBABILITY, MAX_PROBABILITY
export LOD_THRESHOLD_SUGGESTIVE, LOD_THRESHOLD_SIGNIFICANT
export GWAS_GENOME_WIDE_ALPHA, GWAS_SUGGESTIVE_ALPHA
export DEFAULT_MAF_THRESHOLD, DEFAULT_MISSING_RATE_THRESHOLD

# 工具函数
export safe_log, safe_exp, logit, expit, logsumexp
export is_positive_definite, make_positive_definite
export validate_genotypes, validate_phenotypes, validate_population
export load_population, load_pedigree, load_genetic_map, save_results

# 表型与分布
export phenotype_decomposition, variance_partition_phenotypes
export continuous_vs_discrete_traits, phenotypic_correlation
export moments_around_mean, raw_moments, skewness, kurtosis
export cumulants_from_moments, multivariate_normal_dens
export mixture_model_likelihood, transformation_normal

# 回归
export bols_with_intercept, simple_regression, weighted_ols
export parent_offspring_regression, partial_correlation_recursion
export path_coeff_matrix, standardized_coefficients, variance_inflation_factor

# 单位点遗传学
export allele_genotype_freqs, allele_freqs_from_matrix
export hwe_expected_freqs, hwe_chi_square_test, hwe_exact_test
export minor_allele_frequency, heterozygosity, observed_heterozygosity
export fisher_decomposition, average_effects_alpha
export breeding_value_partition, dominance_deviation_calc
export variance_partition_single_locus

# 双位点与LD
export two_locus_epistasis, cockerham_orthogonal_partition
export linkage_disequilibrium_metrics, ld_from_genotypes, ld_matrix

# 环境与GxE
export reaction_norm_slope, gxe_anova

# 亲缘关系
export recursive_kinship, additive_relationship_matrix
export grm_vanraden_yang, realized_relatedness

# 混合模型
export general_linear_model, mixed_model_equations
export animal_model_scan, reml_ai_algorithm

# 品系杂交
export generation_mean_analysis, castle_wright_n_factors
export inbreeding_depression_linear, estimate_lethal_equivalents

# 尺度
export mather_jinks_scaling_test, threshold_model_liability, liability_heritability

# 连锁分析
export haldane_mapping, inverse_haldane, kosambi_mapping, inverse_kosambi
export estimate_recombination_mle

# QTL作图
export haley_knott_scan, interval_mapping_scan, cim_scan
export haseman_elston_regression, outbred_variance_comp_qtl

# GWAS
export gwas_single_marker_scan, gwas_mlm_fit, genomic_control_lambda
export manhattan_plot, qq_plot
export structured_association_mapping, bonferroni_threshold, fdr_correction, ld_clumping

# 基因组预测
export gblup_prediction_solver, rr_blup_marker_effects
export bayesian_abc_models, cross_validation_scheme

# 连锁图谱构建
export construct_genetic_map, two_point_linkage_matrix
export linkage_group_assignment, order_markers_within_group

# 数学工具
export kronecker_product, hadamard_product, generalized_inverse
export vech_operator, duplication_matrix, commutation_matrix
export moment_generating_function, conditional_distribution_mvn
export newton_raphson_update, em_algorithm_template, fisher_scoring_iteration
export beta_function, hermite_polynomial, laguerre_polynomial, bessel_i0

# 类型辅助方法
export is_homozygous, is_heterozygous, genotype_to_numeric
export n_alleles, is_biallelic, n_markers, n_chromosomes
export marker_distance, total_length
export n_observed, trait_mean, trait_var
export n_individuals, n_traits, get_phenotypes, get_genotypes
export is_in_hwe, max_lod, significant_peaks
export dominance_ratio

end # module

