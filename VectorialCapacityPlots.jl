include("ParameterFitting.jl")
include("VCCalculations.jl")

using .ParameterFitting, .VCCalculations

using Plots, LaTeXStrings, Random, Statistics

function rr(i)
    rng = MersenneTwister(i);
    return rng
end

ns = 10000 # number of simulations
CL = 0.95 # confidence for uncertainty from simulations
q = (1-CL)/2 # used for lower and upper (1-q) confidence level

# Useful parameters
α = 0.25; # bite rate
σ = 0.09698706953718507; # incubation rate
λ = 31*σ; # used in Erlang distribution
bites_k = collect(0:50); # bites in step 2
a₀_survEIP = collect(0:64); # ages in step 1
a₁_bites = collect(5:35); # ages in step 2
a₀_bites = collect(0:35); # ages in step 3

μs_AI_control = rand(rr(123), d_AI_control, ns);
μs_AI_treated = rand(rr(123), d_AI_treated, ns);

μs_L_control = exp.(rand(rr(123), d_L_control, ns));
μs_L_treated = exp.(rand(rr(123), d_L_treated, ns));

μs_G_control = rand(rr(123), d_G_control, ns);
μs_G_treated = rand(rr(123), d_G_treated, ns);


# Step 1: What is the probability of a mosquito surviving the EIP given it takes an infectious blood-meal at age a₀?

surviving_EIP_control_ageind = zeros(length(a₀_survEIP), ns)
surviving_EIP_control_logistic = zeros(length(a₀_survEIP), ns)
surviving_EIP_control_gompertz = zeros(length(a₀_survEIP), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    surviving_EIP_control_ageind[:, i] = surviving_EIP(μ, a₀_survEIP, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_control[:, i]
    surviving_EIP_control_logistic[:, i] = surviving_EIP(μ, a₀_survEIP, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_control[:, i]
    surviving_EIP_control_gompertz[:, i] = surviving_EIP(μ, a₀_survEIP, "Gompertz", "Erlang")
end

surviving_EIP_treated_ageind = zeros(length(a₀_survEIP), ns)
surviving_EIP_treated_logistic = zeros(length(a₀_survEIP), ns)
surviving_EIP_treated_gompertz = zeros(length(a₀_survEIP), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    surviving_EIP_treated_ageind[:, i] = surviving_EIP(μ, a₀_survEIP, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_treated[:, i]
    surviving_EIP_treated_logistic[:, i] = surviving_EIP(μ, a₀_survEIP, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_treated[:, i]
    surviving_EIP_treated_gompertz[:, i] = surviving_EIP(μ, a₀_survEIP, "Gompertz", "Erlang")
end

# lower and upper bounds
surviving_EIP_control_ageind_l = [quantile(surviving_EIP_control_ageind[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_control_ageind_u = [quantile(surviving_EIP_control_ageind[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_treated_ageind_l = [quantile(surviving_EIP_treated_ageind[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_treated_ageind_u = [quantile(surviving_EIP_treated_ageind[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_control_logistic_l = [quantile(surviving_EIP_control_logistic[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_control_logistic_u = [quantile(surviving_EIP_control_logistic[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_treated_logistic_l = [quantile(surviving_EIP_treated_logistic[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_treated_logistic_u = [quantile(surviving_EIP_treated_logistic[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_control_gompertz_l = [quantile(surviving_EIP_control_gompertz[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_control_gompertz_u = [quantile(surviving_EIP_control_gompertz[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_treated_gompertz_l = [quantile(surviving_EIP_treated_gompertz[i, :], q) for i in 1:length(a₀_survEIP)]
surviving_EIP_treated_gompertz_u = [quantile(surviving_EIP_treated_gompertz[i, :], 1-q) for i in 1:length(a₀_survEIP)]

surviving_EIP_control_AI = surviving_EIP(μ_ageind_control, a₀_survEIP, "Age-Independent", "Erlang")
surviving_EIP_treated_AI = surviving_EIP(μ_ageind_treated, a₀_survEIP, "Age-Independent", "Erlang")
surviving_EIP_control_L = surviving_EIP(μ_logistic_control, a₀_survEIP, "Logistic", "Erlang")
surviving_EIP_treated_L = surviving_EIP(μ_logistic_treated, a₀_survEIP, "Logistic", "Erlang")
surviving_EIP_control_G = surviving_EIP(μ_gompertz_control, a₀_survEIP, "Gompertz", "Erlang")
surviving_EIP_treated_G = surviving_EIP(μ_gompertz_treated, a₀_survEIP, "Gompertz", "Erlang")


# Step 2: How many bites will the mosquito take if it has survived the EIP?

mean_bites_control_AI = zeros(length(a₁_bites), ns)
mean_bites_control_L = zeros(length(a₁_bites), ns)
mean_bites_control_G = zeros(length(a₁_bites), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    mean_bites_control_AI[:, i] = no_bites(μ, a₁_bites, "Age-Independent").means
end
for i in 1:ns
    μ = μs_L_control[:, i]
    mean_bites_control_L[:, i] = no_bites(μ, a₁_bites, "Logistic").means
end
for i in 1:ns
    μ = μs_G_control[:, i]
    mean_bites_control_G[:, i] = no_bites(μ, a₁_bites, "Gompertz").means
end

bites_PMF_control_AI = zeros(length(bites_k), ns)
bites_PMF_control_L = zeros(length(bites_k), ns)
bites_PMF_control_G = zeros(length(bites_k), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    bites_PMF_control_AI[:, i] = no_bites(μ, a₁_bites[11], "Age-Independent").bites
end
for i in 1:ns
    μ = μs_L_control[:, i]
    bites_PMF_control_L[:, i] = no_bites(μ, a₁_bites[11], "Logistic").bites
end
for i in 1:ns
    μ = μs_G_control[:, i]
    bites_PMF_control_G[:, i] = no_bites(μ, a₁_bites[11], "Gompertz").bites
end

# lower and upper bounds
mean_bites_control_AI_l = [quantile(mean_bites_control_AI[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_control_AI_u = [quantile(mean_bites_control_AI[i, :], 1-q) for i in 1:length(a₁_bites)]

mean_bites_control_L_l = [quantile(mean_bites_control_L[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_control_L_u = [quantile(mean_bites_control_L[i, :], 1-q) for i in 1:length(a₁_bites)]

mean_bites_control_G_l = [quantile(mean_bites_control_G[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_control_G_u = [quantile(mean_bites_control_G[i, :], 1-q) for i in 1:length(a₁_bites)]

bites_PMF_control_AI_l = [quantile(bites_PMF_control_AI[i, :], q) for i in 1:length(bites_k)]
bites_PMF_control_AI_u = [quantile(bites_PMF_control_AI[i, :], 1-q) for i in 1:length(bites_k)]

bites_PMF_control_L_l = [quantile(bites_PMF_control_L[i, :], q) for i in 1:length(bites_k)]
bites_PMF_control_L_u = [quantile(bites_PMF_control_L[i, :], 1-q) for i in 1:length(bites_k)]

bites_PMF_control_G_l = [quantile(bites_PMF_control_G[i, :], q) for i in 1:length(bites_k)]
bites_PMF_control_G_u = [quantile(bites_PMF_control_G[i, :], 1-q) for i in 1:length(bites_k)]

mean_bites_treated_AI = zeros(length(a₁_bites), ns)
mean_bites_treated_L = zeros(length(a₁_bites), ns)
mean_bites_treated_G = zeros(length(a₁_bites), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    mean_bites_treated_AI[:, i] = no_bites(μ, a₁_bites, "Age-Independent").means
end
for i in 1:ns
    μ = μs_L_treated[:, i]
    mean_bites_treated_L[:, i] = no_bites(μ, a₁_bites, "Logistic").means
end
for i in 1:ns
    μ = μs_G_treated[:, i]
    mean_bites_treated_G[:, i] = no_bites(μ, a₁_bites, "Gompertz").means
end

bites_PMF_treated_AI = zeros(length(bites_k), ns)
bites_PMF_treated_L = zeros(length(bites_k), ns)
bites_PMF_treated_G = zeros(length(bites_k), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    bites_PMF_treated_AI[:, i] = no_bites(μ, a₁_bites[11], "Age-Independent").bites
end
for i in 1:ns
    μ = μs_L_treated[:, i]
    bites_PMF_treated_L[:, i] = no_bites(μ, a₁_bites[11], "Logistic").bites
end
for i in 1:ns
    μ = μs_G_treated[:, i]
    bites_PMF_treated_G[:, i] = no_bites(μ, a₁_bites[11], "Gompertz").bites
end

# lower and upper bounds
mean_bites_treated_AI_l = [quantile(mean_bites_treated_AI[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_treated_AI_u = [quantile(mean_bites_treated_AI[i, :], 1-q) for i in 1:length(a₁_bites)]

mean_bites_treated_L_l = [quantile(mean_bites_treated_L[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_treated_L_u = [quantile(mean_bites_treated_L[i, :], 1-q) for i in 1:length(a₁_bites)]

mean_bites_treated_G_l = [quantile(mean_bites_treated_G[i, :], q) for i in 1:length(a₁_bites)]
mean_bites_treated_G_u = [quantile(mean_bites_treated_G[i, :], 1-q) for i in 1:length(a₁_bites)]

bites_PMF_treated_AI_l = [quantile(bites_PMF_treated_AI[i, :], q) for i in 1:length(bites_k)]
bites_PMF_treated_AI_u = [quantile(bites_PMF_treated_AI[i, :], 1-q) for i in 1:length(bites_k)]

bites_PMF_treated_L_l = [quantile(bites_PMF_treated_L[i, :], q) for i in 1:length(bites_k)]
bites_PMF_treated_L_u = [quantile(bites_PMF_treated_L[i, :], 1-q) for i in 1:length(bites_k)]

bites_PMF_treated_G_l = [quantile(bites_PMF_treated_G[i, :], q) for i in 1:length(bites_k)]
bites_PMF_treated_G_u = [quantile(bites_PMF_treated_G[i, :], 1-q) for i in 1:length(bites_k)]

bites_surv_EIP_control_ageind = no_bites(μ_ageind_control, a₁_bites, "Age-Independent").bites;
means_bites_control_ageind = no_bites(μ_ageind_control, a₁_bites, "Age-Independent").means;

bites_surv_EIP_control_logistic = no_bites(μ_logistic_control, a₁_bites, "Logistic").bites;
means_bites_control_logistic = no_bites(μ_logistic_control, a₁_bites, "Logistic").means;

bites_surv_EIP_control_gompertz = no_bites(μ_gompertz_control, a₁_bites, "Gompertz").bites;
means_bites_control_gompertz = no_bites(μ_gompertz_control, a₁_bites, "Gompertz").means;

bites_surv_EIP_treated_ageind = no_bites(μ_ageind_treated, a₁_bites, "Age-Independent").bites;
means_bites_treated_ageind = no_bites(μ_ageind_treated, a₁_bites, "Age-Independent").means;

bites_surv_EIP_treated_logistic = no_bites(μ_logistic_treated, a₁_bites, "Logistic").bites;
means_bites_treated_logistic = no_bites(μ_logistic_treated, a₁_bites, "Logistic").means;

bites_surv_EIP_treated_gompertz = no_bites(μ_gompertz_treated, a₁_bites, "Gompertz").bites;
means_bites_treated_gompertz = no_bites(μ_gompertz_treated, a₁_bites, "Gompertz").means;

mean_error_control_AI = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_control_AI[i] = (means_bites_control_ageind[i] - mean_bites_control_AI_l[i], -(means_bites_control_ageind[i] - mean_bites_control_AI_u[i]))
end

mean_error_control_L = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_control_L[i] = (means_bites_control_logistic[i] - mean_bites_control_L_l[i], -(means_bites_control_logistic[i] - mean_bites_control_L_u[i]))
end

mean_error_control_G = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_control_G[i] = (means_bites_control_gompertz[i] - mean_bites_control_G_l[i], -(means_bites_control_gompertz[i] - mean_bites_control_G_u[i]))
end

mean_error_treated_AI = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_treated_AI[i] = (means_bites_treated_ageind[i] - mean_bites_treated_AI_l[i], -(means_bites_treated_ageind[i] - mean_bites_treated_AI_u[i]))
end

mean_error_treated_L = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_treated_L[i] = (means_bites_treated_logistic[i] - mean_bites_treated_L_l[i], -(means_bites_treated_logistic[i] - mean_bites_treated_L_u[i]))
end

mean_error_treated_G = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    mean_error_treated_G[i] = (means_bites_treated_gompertz[i] - mean_bites_treated_G_l[i], -(means_bites_treated_gompertz[i] - mean_bites_treated_G_u[i]))
end

PMF_error_control_AI = Array{Tuple{Float64, Float64}}(undef, length(bites_k))
for i in 1:length(bites_k)
    PMF_error_control_AI[i] = (bites_surv_EIP_control_ageind[i, 11] - bites_PMF_control_AI_l[i], -(bites_surv_EIP_control_ageind[i, 11] - bites_PMF_control_AI_u[i]))
end

PMF_error_control_L = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    PMF_error_control_L[i] = (bites_surv_EIP_control_logistic[i, 11] - bites_PMF_control_L_l[i], -(bites_surv_EIP_control_logistic[i, 11] - bites_PMF_control_L_u[i]))
end

PMF_error_control_G = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    PMF_error_control_G[i] = (bites_surv_EIP_control_gompertz[i, 11] - bites_PMF_control_G_l[i], -(bites_surv_EIP_control_gompertz[i, 11] - bites_PMF_control_G_u[i]))
end

PMF_error_treated_AI = Array{Tuple{Float64, Float64}}(undef, length(bites_k))
for i in 1:length(bites_k)
    PMF_error_treated_AI[i] = (bites_surv_EIP_treated_ageind[i, 11] - bites_PMF_treated_AI_l[i], -(bites_surv_EIP_treated_ageind[i, 11] - bites_PMF_treated_AI_u[i]))
end

PMF_error_treated_L = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    PMF_error_treated_L[i] = (bites_surv_EIP_treated_logistic[i, 11] - bites_PMF_treated_L_l[i], -(bites_surv_EIP_treated_logistic[i, 11] - bites_PMF_treated_L_u[i]))
end

PMF_error_treated_G = Array{Tuple{Float64, Float64}}(undef, length(a₁_bites))
for i in 1:length(a₁_bites)
    PMF_error_treated_G[i] = (bites_surv_EIP_treated_gompertz[i, 11] - bites_PMF_treated_G_l[i], -(bites_surv_EIP_treated_gompertz[i, 11] - bites_PMF_treated_G_u[i]))
end


# Step 3: What is the expected number of bites a mosquito takes in its lifetime if it has taken an infectious blood-meal at age a₀?

bites_age_control_AI = zeros(length(a₀_bites), ns)
bites_age_control_L = zeros(length(a₀_bites), ns)
bites_age_control_G = zeros(length(a₀_bites), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    bites_age_control_AI[:, i] = bites_age(μ, a₀_bites, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_control[:, i]
    bites_age_control_L[:, i] = bites_age(μ, a₀_bites, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_control[:, i]
    bites_age_control_G[:, i] = bites_age(μ, a₀_bites, "Gompertz", "Erlang")
end

bites_age_treated_AI = zeros(length(a₀_bites), ns)
bites_age_treated_L = zeros(length(a₀_bites), ns)
bites_age_treated_G = zeros(length(a₀_bites), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    bites_age_treated_AI[:, i] = bites_age(μ, a₀_bites, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_treated[:, i]
    bites_age_treated_L[:, i] = bites_age(μ, a₀_bites, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_treated[:, i]
    bites_age_treated_G[:, i] = bites_age(μ, a₀_bites, "Gompertz", "Erlang")
end

bites_age_control_AI_l = [quantile(bites_age_control_AI[i, :], q) for i in 1:length(a₀_bites)]
bites_age_control_L_l = [quantile(bites_age_control_L[i, :], q) for i in 1:length(a₀_bites)]
bites_age_control_G_l = [quantile(bites_age_control_G[i, :], q) for i in 1:length(a₀_bites)]

bites_age_control_AI_u = [quantile(bites_age_control_AI[i, :],1-q) for i in 1:length(a₀_bites)]
bites_age_control_L_u = [quantile(bites_age_control_L[i, :],1-q) for i in 1:length(a₀_bites)]
bites_age_control_G_u = [quantile(bites_age_control_G[i, :],1-q) for i in 1:length(a₀_bites)]

bites_age_treated_AI_l = [quantile(bites_age_treated_AI[i, :], q) for i in 1:length(a₀_bites)]
bites_age_treated_L_l = [quantile(bites_age_treated_L[i, :], q) for i in 1:length(a₀_bites)]
bites_age_treated_G_l = [quantile(bites_age_treated_G[i, :], q) for i in 1:length(a₀_bites)]

bites_age_treated_AI_u = [quantile(bites_age_treated_AI[i, :],1-q) for i in 1:length(a₀_bites)]
bites_age_treated_L_u = [quantile(bites_age_treated_L[i, :],1-q) for i in 1:length(a₀_bites)]
bites_age_treated_G_u = [quantile(bites_age_treated_G[i, :],1-q) for i in 1:length(a₀_bites)]

bites_age_control_ageind = bites_age(μ_ageind_control, a₀_bites, "Age-Independent", "Erlang");
bites_age_control_logistic = bites_age(μ_logistic_control, a₀_bites, "Logistic", "Erlang");
bites_age_control_gompertz = bites_age(μ_gompertz_control, a₀_bites, "Gompertz", "Erlang");

bites_age_treated_ageind = bites_age(μ_ageind_treated, a₀_bites, "Age-Independent", "Erlang");
bites_age_treated_logistic = bites_age(μ_logistic_treated, a₀_bites, "Logistic", "Erlang");
bites_age_treated_gompertz = bites_age(μ_gompertz_treated, a₀_bites, "Gompertz", "Erlang");


# Step 4: What is the expected number of (infectious) bites a mosquito takes in its lifetime?

expected_bites_control_AI = zeros(ns)
expected_bites_control_L = zeros(ns)
expected_bites_control_G = zeros(ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    expected_bites_control_AI[i] = expected_bites(μ, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_control[:, i]
    expected_bites_control_L[i] = expected_bites(μ, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_control[:, i]
    expected_bites_control_G[i] = expected_bites(μ, "Gompertz", "Erlang")
end

expected_bites_treated_AI = zeros(ns)
expected_bites_treated_L = zeros(ns)
expected_bites_treated_G = zeros(ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    expected_bites_treated_AI[i] = expected_bites(μ, "Age-Independent", "Erlang")
end
for i in 1:ns
    μ = μs_L_treated[:, i]
    expected_bites_treated_L[i] = expected_bites(μ, "Logistic", "Erlang")
end
for i in 1:ns
    μ = μs_G_treated[:, i]
    expected_bites_treated_G[i] = expected_bites(μ, "Gompertz", "Erlang")
end

expected_bites_control_AI_l = quantile(expected_bites_control_AI, q)
expected_bites_control_L_l = quantile(expected_bites_control_L, q)
expected_bites_control_G_l = quantile(expected_bites_control_G, q)
expected_bites_control_AI_u = quantile(expected_bites_control_AI, 1-q)
expected_bites_control_L_u = quantile(expected_bites_control_L, 1-q)
expected_bites_control_G_u = quantile(expected_bites_control_G, 1-q)

expected_bites_treated_AI_l = quantile(expected_bites_treated_AI, q)
expected_bites_treated_L_l = quantile(expected_bites_treated_L, q)
expected_bites_treated_G_l = quantile(expected_bites_treated_G, q)
expected_bites_treated_AI_u = quantile(expected_bites_treated_AI, 1-q)
expected_bites_treated_L_u = quantile(expected_bites_treated_L, 1-q)
expected_bites_treated_G_u = quantile(expected_bites_treated_G, 1-q)


expected_bites_ageind_control = expected_bites(μ_ageind_control, "Age-Independent", "Erlang")
expected_bites_logistic_control = expected_bites(μ_logistic_control, "Logistic", "Erlang")
expected_bites_gompertz_control = expected_bites(μ_gompertz_control, "Gompertz", "Erlang")

expected_bites_ageind_treated = expected_bites(μ_ageind_treated, "Age-Independent", "Erlang")
expected_bites_logistic_treated = expected_bites(μ_logistic_treated, "Logistic", "Erlang")
expected_bites_gompertz_treated = expected_bites(μ_gompertz_treated, "Gompertz", "Erlang")


# Decrease in VC between control and treatment (as explained in manuscript)
ageinderl = TreatmentDecrease(expected_bites_ageind_control, expected_bites_ageind_treated)

logisticerl = TreatmentDecrease(expected_bites_logistic_control, expected_bites_logistic_treated)

gompertzerl = TreatmentDecrease(expected_bites_gompertz_control, expected_bites_gompertz_treated)


# Plots for paper

probsurvEIPcontrol = plot(a₀_survEIP, surviving_EIP_control_ageind_l, fillrange=surviving_EIP_control_ageind_u, l=0, fillalpha=0.1, c=:purple4, label=:false);
plot!(a₀_survEIP, surviving_EIP_control_AI, legend=:topright, line=(3, :dot), xlabel=L"a_0~\textrm{(age infectious blood-meal taken)}", label=:false, color=:purple4, ylabel=L"P(\textrm{surviving EIP}~|~\textrm{blood-meal at}~a_0)", legendtitle="Control", foreground_color_legend=nothing, legendtitlefontsize=15, show=true);
plot!(a₀_survEIP, surviving_EIP_control_logistic_l, fillrange=surviving_EIP_control_logistic_u, l=0, fillalpha=0.1, c=:teal, label=:false);
plot!(a₀_survEIP, surviving_EIP_control_L, line=(3, :dash), label=:false, color=:teal);
plot!(a₀_survEIP, surviving_EIP_control_gompertz_l, fillrange=surviving_EIP_control_gompertz_u, l=0, fillalpha=0.1, c=:orange, label=:false);
plot!(a₀_survEIP, surviving_EIP_control_G, line=(3), label=:false, color=:orange);

probsurvEIPtreated = plot(a₀_survEIP, surviving_EIP_treated_ageind_l, fillrange=surviving_EIP_treated_ageind_u, l=0, fillalpha=0.1, c=:purple4, label=:false);
plot!(a₀_survEIP, surviving_EIP_treated_AI, legend=:topright, line=(3, :dot), xlabel=L"a_0~\textrm{(age infectious blood-meal taken)}", label=:false, color=:purple4, ylabel=L"P(\textrm{surviving EIP}~|~\textrm{blood-meal at}~a_0)", legendtitle="Treated", foreground_color_legend=nothing, legendtitlefontsize=15, show=true);
plot!(a₀_survEIP, surviving_EIP_treated_logistic_l, fillrange=surviving_EIP_treated_logistic_u, l=0, fillalpha=0.1, c=:teal, label=:false);
plot!(a₀_survEIP, surviving_EIP_treated_L, line=(3, :dash), label=:false, color=:teal);
plot!(a₀_survEIP, surviving_EIP_treated_gompertz_l, fillrange=surviving_EIP_treated_gompertz_u, l=0, fillalpha=0.1, c=:orange, label=:false);
plot!(a₀_survEIP, surviving_EIP_treated_G, line=(3), label=:false, color=:orange);

probsurvEIP = [probsurvEIPcontrol, probsurvEIPtreated]

legendplot = plot([], label="Gompertz", color=:orange, grid=false, showaxis=false, line=3, legend=:outertop, legendfontsize=10, ticks=false);
plot!([], label="Age-Independent", color=:purple4, line=(3, :dot), ticks=false);
plot!([], label="Logistic", color=:teal, line=(3, :dash), foreground_color_legend=nothing, ticks=false);

l = @layout[grid(1, 2)
            b{0.1h}]

step1 = plot(probsurvEIP..., legendplot, layout=l, ylims=(0,1.0), xlims=(0, 64), xticks=(0:8:64), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10);
display(step1)

# Step 2
mean_ageind_control = round(means_bites_control_ageind[11], sigdigits=3);
mean_ageind_control_ci = round.((mean_bites_control_AI_l[11], mean_bites_control_AI_u[11]), sigdigits=3);
mean_ageind_treated = round(means_bites_treated_ageind[11], sigdigits=3);
mean_ageind_treated_ci = round.((mean_bites_treated_AI_l[11], mean_bites_treated_AI_u[11]), sigdigits=3);
mean_logistic_control = round(means_bites_control_logistic[11], sigdigits=3);
mean_logistic_control_ci = round.((mean_bites_control_L_l[11], mean_bites_control_L_u[11]), sigdigits=3);
mean_logistic_treated = round(means_bites_treated_logistic[11], sigdigits=3);
mean_logistic_treated_ci = round.((mean_bites_treated_L_l[11], mean_bites_control_L_u[11]), sigdigits=3);
mean_gompertz_control = round(means_bites_control_gompertz[11], sigdigits=3);
mean_gompertz_control_ci = round.((mean_bites_control_G_l[11], mean_bites_control_G_u[11]), sigdigits=3);
mean_gompertz_treated = round(means_bites_treated_gompertz[11], sigdigits=3);
mean_gompertz_treated_ci = round.((mean_bites_treated_G_l[11], mean_bites_treated_G_u[11]), sigdigits=3);

pmf_ageind = scatter(bites_k, bites_surv_EIP_control_ageind[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), yerr=PMF_error_control_AI, label="control - $mean_ageind_control $mean_ageind_control_ci");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_ageind[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, yerr=PMF_error_treated_AI, label="treated - $mean_ageind_treated $mean_ageind_treated_ci", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=j~|\textrm{~exits EIP at~}a_1 = 15)");

pmf_logistic = scatter(bites_k, bites_surv_EIP_control_logistic[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), yerr=PMF_error_control_L, label="control - $mean_logistic_control $mean_logistic_control_ci");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_logistic[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, yerr=PMF_error_treated_L, label="treated - $mean_logistic_treated $mean_logistic_treated_ci", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=j~|\textrm{~exits EIP at~}a_1 = 15)");

pmf_gompertz = scatter(bites_k, bites_surv_EIP_control_gompertz[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), yerr=PMF_error_control_G, label="control - $mean_gompertz_control $mean_gompertz_control_ci");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_gompertz[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, yerr=PMF_error_treated_G, label="treated - $mean_gompertz_treated $mean_gompertz_treated_ci", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=j~|\textrm{~exits EIP at~}a_1 = 15)", xlabel=L"\textrm{Number of bites}~(j)");

pmfs = [pmf_ageind, pmf_logistic, pmf_gompertz];

pmfs_panel = plot(pmfs..., layout=grid(3,1), fontfamily="Times", size=(600, 1000), show=true, tickfontsize=10, xguidefontsize=12);
display(pmfs_panel)


# Heatmap panel
heatmap_ageind_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), ylabel="Control", legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_ageind, color=:gray, marker=:square, label="", alpha=0.75, yerr=mean_error_control_AI, msc=:gray);

heatmap_ageind_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), ylabel="Treated", legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_ageind, color=:gray, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}", yerr=mean_error_treated_AI, msc=:gray);

heatmap_logistic_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_logistic, color=:gray, marker=:square, label="", alpha=0.75, yerr=mean_error_control_L, msc=:gray);

heatmap_logistic_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_logistic, color=:gray, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}", yerr=mean_error_treated_L, msc=:gray);

heatmap_gompertz_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_gompertz, color=:gray, marker=:square, label="", alpha=0.75, yerr=mean_error_control_G, msc=:gray);

heatmap_gompertz_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_gompertz, color=:gray, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}", yerr=mean_error_treated_G, msc=:gray);

heatmaps = [heatmap_ageind_control, heatmap_logistic_control, heatmap_gompertz_control, heatmap_ageind_treated, heatmap_logistic_treated, heatmap_gompertz_treated]

legendplot = scatter([], showaxis = false, grid = false, color=:gray, alpha=0.75, marker=:square, label="mean number of bites", msc=:gray, foreground_color_legend=nothing, legend=:outertop, legendfontsize=10, ticks=false);

l = @layout[grid(2,3)
            e{0.05h}]

heatmap_panel = plot(heatmaps..., legendplot, layout=l, fontfamily="Times", size=(1900, 1300), show=true, tickfontsize=15, xguidefontsize=20, yticks=0:3:12);
display(heatmap_panel)


# Step 3 plots

expected_bites_control = plot(a₀_bites, bites_age_control_AI_l, fillrange=bites_age_control_AI_u, fillalpha=0.1, line=0, label=:false, c=:purple4);
plot!(a₀_bites, bites_age_control_ageind, line=(3, :dot), c=:purple4, legend=:false, xlabel=L"a_0~\textrm{(age infectious blood-meal taken in days)}", ylabel=L"E(z~|~\textrm{blood-meal at}~a_0)");
plot!(a₀_bites, bites_age_control_L_l, fillrange=bites_age_control_L_u, fillalpha=0.1, line=0, c=:teal);
plot!(a₀_bites, bites_age_control_logistic, line=(3, :dash), c=:teal);
plot!(a₀_bites, bites_age_control_G_l, fillrange=bites_age_control_G_u, fillalpha=0.1, l=0, c=:orange);
plot!(a₀_bites, bites_age_control_gompertz, line=3, c=:orange);

expected_bites_treated = plot(a₀_bites, bites_age_treated_AI_l, fillrange=bites_age_treated_AI_u, fillalpha=0.1, line=0, label=:false, c=:purple4);
plot!(a₀_bites, bites_age_treated_ageind, line=(3, :dot), c=:purple4, legend=:false, xlabel=L"a_0~\textrm{(age infectious blood-meal taken in days)}", ylabel=L"E(z~|~\textrm{blood-meal at}~a_0)");
plot!(a₀_bites, bites_age_treated_L_l, fillrange=bites_age_treated_L_u, fillalpha=0.1, line=0, c=:teal);
plot!(a₀_bites, bites_age_treated_logistic, line=(3, :dash), c=:teal);
plot!(a₀_bites, bites_age_treated_G_l, fillrange=bites_age_treated_G_u, fillalpha=0.1, l=0, c=:orange);
plot!(a₀_bites, bites_age_treated_gompertz, line=3, c=:orange);

expected_bites_plots = [expected_bites_control, expected_bites_treated]

legendplot = plot([], label="Gompertz", color=:orange, grid=false, showaxis=false, line=3, legend=:outertop, legendfontsize=10);
plot!([], label="Age-Independent", color=:purple4, line=(3, :dot));
plot!([], label="Logistic", color=:teal, line=(3, :dash), foreground_color_legend=nothing, ticks=false);

l = @layout[grid(1, 2)
            b{0.1h}]
step3 = plot(expected_bites_plots..., legendplot, layout=l, ylims=(0, 6), xticks=(0:5:35), xlims=(0, 35), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10, show=true);
display(step3)


# Step 4 plots

bites_error_control_AI = [(expected_bites_ageind_control - expected_bites_control_AI_l, -(expected_bites_ageind_control - expected_bites_control_AI_u))]
bites_error_control_L =[(expected_bites_logistic_control - expected_bites_control_L_l, -(expected_bites_logistic_control - expected_bites_control_L_u))]
bites_error_control_G =[(expected_bites_gompertz_control - expected_bites_control_G_l, -(expected_bites_gompertz_control - expected_bites_control_G_u))]
bites_error_treated_AI = [(expected_bites_ageind_treated - expected_bites_treated_AI_l, -(expected_bites_ageind_treated - expected_bites_treated_AI_u))]
bites_error_treated_L =[(expected_bites_logistic_treated - expected_bites_treated_L_l, -(expected_bites_logistic_treated - expected_bites_treated_L_u))]
bites_error_treated_G =[(expected_bites_gompertz_treated - expected_bites_treated_G_l, -(expected_bites_gompertz_treated - expected_bites_treated_G_u))]

bites_control = scatter([1], [expected_bites_ageind_control], xlims=(0,4), ylims=(0,6), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", show=true, ms=5, m=:rect, label="", c=:purple4, title="Control", size=(400,380), ylabel=("Expected number of bites"), yerr=bites_error_control_AI, msc=:purple4);
scatter!([2], [expected_bites_logistic_control], ms=5, m=:rect, label="", c=:teal, yerr=bites_error_control_L, msc=:teal);
scatter!([3], [expected_bites_gompertz_control], ms=5, m=:rect, label="", c=:orange, tickdirection=:out, yerr=bites_error_control_G, msc=:orange);

bites_treated = scatter([1], [expected_bites_ageind_treated], xlims=(0,4), ylims=(0,6), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", show=true, ms=5, m=:rect, label="", c=:purple4, title="Treated", size=(400,380), yerr=bites_error_treated_AI, msc=:purple4);
scatter!([2], [expected_bites_logistic_treated], ms=5, m=:rect, label="", c=:teal, yerr=bites_error_treated_L, msc=:teal);
scatter!([3], [expected_bites_gompertz_treated], ms=5, m=:rect, label="", c=:orange, tickdirection=:out, yerr=bites_error_treated_G, msc=:orange);

bites_plots = [bites_control, bites_treated]

step4 = plot(bites_plots..., layout=grid(1,2), ylims=(0, 6), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10, show=true);
display(step4)


# Violin plot (final figure)
using StatsPlots

VC_AI = TreatmentDecrease(expected_bites_control_AI, expected_bites_treated_AI)
VC_L = TreatmentDecrease(expected_bites_control_L, expected_bites_treated_L)
VC_G = TreatmentDecrease(expected_bites_control_G, expected_bites_treated_G)

VC_AI_l = quantile(VC_AI, q)
VC_AI_u = quantile(VC_AI, 1-q)
VC_L_l = quantile(VC_L, q)
VC_L_u = quantile(VC_L, 1-q)
VC_G_l = quantile(VC_G, q)
VC_G_u = quantile(VC_G, 1-q)

violinplot = violin(repeat([1], 10000), VC_AI, alpha=0.3, legend=:none, tickdirection=:out, c=:purple4, ylims=(0, 100), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", ylabel="Percentage relative to control");
scatter!([1], [ageinderl], c=:purple4, ms=3, m=:rect, yerr=[(ageinderl-VC_AI_l, -(ageinderl-VC_AI_u))], msc=:purple4);
violin!(repeat([2], 10000), VC_L, alpha=0.3, c=:teal);
scatter!([2], [logisticerl], c=:teal, ms=3, m=:rect, yerr=[(logisticerl-VC_L_l, -(logisticerl-VC_L_u))], msc=:teal);
violin!(repeat([3], 10000), VC_G, alpha=0.3, c=:orange);
scatter!([3], [gompertzerl], c=:orange, ms=3, m=:rect, yerr=[(gompertzerl-VC_G_l, -(gompertzerl-VC_G_u))], msc=:orange);
display(violinplot)
