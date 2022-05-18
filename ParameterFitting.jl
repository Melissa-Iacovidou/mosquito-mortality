module ParameterFitting

export ages_control, ages_treated, Mortality, Survival1,  μ_ageind_control, μ_ageind_control_error, μ_ageind_treated, μ_ageind_treated_error, μ_logistic_control, μ_logistic_control_error, μ_logistic_treated, μ_logistic_treated_error, μ_gompertz_control, μ_gompertz_control_error, μ_gompertz_treated, μ_gompertz_treated_error, KM_control, KM_treated, CIU_C, CIU_T, CIL_C, CIL_T, d_AI_control, d_AI_treated, d_G_control, d_G_treated, d_L_control, d_L_treated

using XLSX, DataFrames, Survival, Optim, Distributions, Measurements, LinearAlgebra

# Edit file paths accordingly

# Load survival data - CONTROL
df_survival_control = DataFrame(XLSX.readtable("../../Supporting_information/S5_Table.xlsx", "SurvivalControl")...);
df_KM_control = DataFrame(XLSX.readtable("../../Supporting_Information/S7_Table.xlsx", "KMControl")...);

# Load survival data - TREATED
df_survival_treated = DataFrame(XLSX.readtable("../../Supporting_information/S6_Table.xlsx", "SurvivalTreated")...);
df_KM_treated = DataFrame(XLSX.readtable("../../Supporting_Information/S8_Table.xlsx", "KMTreated")...);

# Kaplan-Meier Estimator (using Survival.jl package)
# Breaking down the upper and lower confidence interval for plotting purposes
status = ones(Int64, 200);

times_control = convert(Array{Int64}, df_KM_control.ti);
times_treated = convert(Array{Int64}, df_KM_treated.ti);

KM_control = Survival.fit(KaplanMeier, times_control, status);
CIL_C = zeros(length(confint(KM_control)))
CIU_C = zeros(length(confint(KM_control)))
for i in 1:length(CIL_C)
    CIL_C[i] = confint(KM_control)[i][1]
    CIU_C[i] = confint(KM_control)[i][2]
end

KM_treated = Survival.fit(KaplanMeier, times_treated, status);
CIL_T = zeros(length(confint(KM_treated)))
CIU_T = zeros(length(confint(KM_treated)))
for i in 1:length(CIL_T)
    CIL_T[i] = confint(KM_treated)[i][1]
    CIU_T[i] = confint(KM_treated)[i][2]
end

# Ages from data
ages_control = collect(skipmissing(df_survival_control.Age));
ages_treated = collect(skipmissing(df_survival_treated.Age));

# Initial parameter values
μ₀_ageind = [0.05];
μ₀_logistic = [-1.5, -2.0, 3.0]; # note transformation in S1_File
μ₀_gompertz = [0.005, 0.1];

# Mortality function
function Mortality(a, μ, case)
    if case == "Age-Independent"
        if length(μ) != 1
            throw(ArgumentError("Length of μ must be = 1 for Age-Independent case"))
        else
            return μ .* ones(length(a))
        end
    elseif case == "Logistic"
        if length(μ) != 3
            throw(ArgumentError("Length of μ must be = 3 for Logistic case"))
        else
            return μ[1] ./ (1 .+ exp.(μ[2] .* (-a .+ μ[3])))
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            throw(ArgumentError("Length of μ must be = 2 for Gompertz case"))
        else
            return μ[1] .* exp.(μ[2] .* a)
        end
    else
        throw(ArgumentError("Case must be either Age-Independent, Logistic, or Gompertz"))
    end
end

# Survival function (named Survival1 dues to conflicts with Survival.jl package)
function Survival1(a, μ, case)
    if case == "Age-Independent"
        if length(μ) != 1
            throw(ArgumentError("Length of μ must be = 1 for Age-Independent case"))
        else
            return exp.(-μ.*a)
        end
    elseif case == "Logistic"
        if length(μ) != 3
            throw(ArgumentError("Length of μ must be = 3 for Logistic case"))
        else
            return exp.(-a.*μ[1] .+ μ[1]/μ[2] .* log.(1 .+ exp.(μ[2] * μ[3])) .- μ[1]/μ[2] .* log.(1 .+ exp.(μ[2].*(μ[3] .- a))))
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            throw(ArgumentError("Length of μ must be = 2 for Gompertz case"))
        else
            return exp.(-μ[1] .* (exp.(μ[2] .* a) .- 1) ./ μ[2])
        end
    else
        throw(ArgumentError("Case must be either Age-Independent, Logistic, or Gompertz"))
    end
end

# Maximum Likelihood Estimation Fitting (S1_File for details)
# Define negative log-likelihood functions for all cases
# Log Transform all parameters in Logistic case
function NegLogLikelihood(μ, case, times; k=200)
    if case == "Age-Independent" && length(μ) != 1 || case == "Logistic" && length(μ) != 3 || case == "Gompertz" && length(μ) != 2
        throw(ArgumentError("Length of μ/params_iv doesn't match case."))
    end

    if case == "Age-Independent"
        return - (k*log(μ[1]) - μ[1]*sum(times))
    elseif case == "Logistic"
        return - (k*μ[1] + k*exp(μ[1] - μ[2]) * log(1.0+exp(exp(μ[2]+μ[3]))) - exp(μ[1])*sum(times) - (exp(μ[1]-μ[2]) + 1.0) * sum(log.(1.0 .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- times)))))

    elseif case == "Gompertz"
        return - (k*log(μ[1]) + μ[2]*sum(times) + k*μ[1]/μ[2] - μ[1]/μ[2]*sum(exp.(μ[2] .* times)))
    else
        throw(ArgumentError("Case not supported"))
    end
end

# Obtaining the MLE (using Optim.jl package and NegLogLikelihood function)
# Note: Logistic case we allow negative values due to the transformation of the variables (see S1_File)
function fit(case, times, params_iv)

    function NLL(μ)
        return NegLogLikelihood(μ, case, times)
    end
    if case == "Logistic"
        mle = optimize(NLL, -Inf, Inf, params_iv, Fminbox(LBFGS()))
    else
        mle = optimize(NLL, 0, Inf, params_iv, Fminbox(LBFGS()))
    end
    return mle
end

# Fitting results
# Age-Independent
fit_control_ageind = fit("Age-Independent", times_control, μ₀_ageind);
fit_treated_ageind = fit("Age-Independent", times_treated, μ₀_ageind);

μ_ageind_control = fit_control_ageind.minimizer;
μ_ageind_treated = fit_treated_ageind.minimizer;

# Logistic
fit_control_logistic = fit("Logistic", times_control, μ₀_logistic);
fit_treated_logistic = fit("Logistic", times_treated, μ₀_logistic);

μ_logistic_control = fit_control_logistic.minimizer;
μ_logistic_treated = fit_treated_logistic.minimizer;

# Gompertz
fit_control_gompertz = fit("Gompertz", times_control, μ₀_gompertz);
fit_treated_gompertz = fit("Gompertz", times_treated, μ₀_gompertz);

μ_gompertz_control = fit_control_gompertz.minimizer;
μ_gompertz_treated = fit_treated_gompertz.minimizer;

# Fisher Information Matrix for multiple parameters (including transformation in Logistic case) - matrix comprises of second derivatives of log-likelihood functions - returns inverse of FIM used in Wald confidence intervals (see S1_File)
function FisherInformationMatrix(case, μ, aᵢ; k=200)
    if case == "Logistic"
        D2LDM12 = k*exp(μ[1]-μ[2])*log(1+exp(exp(μ[2]+μ[3]))) - exp(μ[1])*sum(aᵢ) - exp(μ[1]-μ[2])*sum(log.(1 .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- aᵢ))))

        D2LDM1DM2 = - k*exp(μ[1]-μ[2])*log(1+exp(exp(μ[2]+μ[3]))) + k*exp(μ[1]+μ[3])/(1+exp(-exp(μ[2]+μ[3]))) + exp(μ[1]-μ[2])*sum(log.(1 .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- aᵢ)))) - exp(μ[1])*sum((exp(μ[3]) .- aᵢ) ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ))))

        D2LDM1DM3 = k*exp(μ[1]+μ[3])/(1+exp(-exp(μ[2]+μ[3]))) - exp(μ[1]+μ[3])*sum(1 ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ))))

        D2LDM2DM1 = D2LDM1DM2

        D2LDM22 = k*exp(μ[1]-μ[2]) * log(1+exp(exp(μ[2]+μ[3]))) - k*exp(μ[1]+μ[3])*exp(-exp(μ[2]+μ[3]))*(1-exp(μ[2]+μ[3])+exp(exp(μ[2]+μ[3])))/((1+exp(-exp(μ[2]+μ[3])))^2) - exp(μ[1]-μ[2]) * sum(log.(1 .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- aᵢ)))) + (exp(μ[1])-exp(μ[2])) * sum((exp(μ[3]) .- aᵢ) ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ)))) - exp(μ[2])*(exp(μ[1]) + exp(μ[2])) * sum((exp(μ[3]) .- aᵢ).^2 .* exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ)) ./ ((1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ))).^2))

        D2LDM2DM3 = k*exp(μ[1]+μ[2] + 2*μ[3])*exp(-exp(μ[2]+μ[3])) /(1+exp(-exp(μ[2]+μ[3])))^2 + exp(μ[1]+μ[3])*sum(1 ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ)))) - (exp(μ[1]+μ[3])+exp(μ[2]+μ[3]))*sum(exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ)) .* (1 .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- aᵢ)) .+ exp(μ[2]+μ[3]) .- aᵢ .* exp(μ[2])) ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ))).^2)

        D2LDM3DM1 = D2LDM1DM3

        D2LDM3DM2 = D2LDM2DM3

        D2LDM32 = k*exp(μ[1]+μ[3])*exp(-exp(μ[2]+μ[3])) * (1+exp(exp(μ[2]+μ[3]))+exp(μ[2]+μ[3])) / (1+exp(-exp(μ[2]+μ[3])))^2 - (exp(μ[1]+μ[3])+exp(μ[2]+μ[3]))*sum(exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ)) .* (1 .+ exp(μ[2]+μ[3]) .+ exp.(exp(μ[2]) .* (exp(μ[3]) .- aᵢ))) ./ (1 .+ exp.(-exp(μ[2]) .* (exp(μ[3]) .- aᵢ))).^2)

        FIM = - [D2LDM12 D2LDM1DM2 D2LDM1DM3
                 D2LDM2DM1 D2LDM22 D2LDM2DM3
                 D2LDM3DM1 D2LDM3DM2 D2LDM32]

        FIM_inv = inv(FIM)
    elseif case == "Gompertz"
        D2LDG12 = - k / μ[1]^2
        D2LDG1DG2 = - k / μ[2]^2 + 1/μ[2]^2 * sum(exp.(μ[2] .* aᵢ)) - 1/μ[2] * sum(aᵢ .* exp.(μ[2] .* aᵢ))
        D2LDG2DG1 = D2LDG1DG2
        D2LDG22 = 2*k*μ[1]/(μ[2]^3) - 2*μ[1]/(μ[2]^3)*sum(exp.(μ[2] .* aᵢ)) + 2*μ[1]/(μ[2]^2) * sum(aᵢ .* exp.(μ[2] .* aᵢ)) - μ[1]/μ[2] * sum(aᵢ.^2 .* exp.(μ[2] .* aᵢ))

        FIM = - [D2LDG12 D2LDG1DG2
                 D2LDG2DG1 D2LDG22]

        FIM_inv = inv(FIM)
    else
        throw(ArgumentError("Case invalid"))
    end
    return FIM_inv
end

z = quantile(Normal(), 0.975) # 95% critical value
k = 200 # number of samples

# Calculating the Wald confidence intervals
# Age-Independent
FIM_AI_control = μ_ageind_control.^2 ./ k  #inverse of negative second derivative of log-likelihood function used for confidence interval calculation (see S1_File)
FIM_AI_treated = μ_ageind_treated.^2 ./ k

μ_ageind_control_error = μ_ageind_control .± z .* sqrt.(FIM_AI_control)
μ_ageind_treated_error = μ_ageind_treated .± z .* sqrt.(FIM_AI_treated)

# Logistic
FIM_L_control = FisherInformationMatrix("Logistic", μ_logistic_control, times_control)
FIM_L_treated = FisherInformationMatrix("Logistic", μ_logistic_treated, times_treated)

μ_logistic_control_error = μ_logistic_control .± z .* sqrt.([FIM_L_control[1,1], FIM_L_control[2,2], FIM_L_control[3,3]])
μ_logistic_treated_error = μ_logistic_treated .± z .* sqrt.([FIM_L_treated[1,1], FIM_L_treated[2,2], FIM_L_treated[3,3]])

# Gompertz
FIM_G_control = FisherInformationMatrix("Gompertz", μ_gompertz_control, times_control)
FIM_G_treated = FisherInformationMatrix("Gompertz", μ_gompertz_treated, times_treated)

μ_gompertz_control_error = μ_gompertz_control .± z .* sqrt.([FIM_G_control[1,1], FIM_G_control[2,2]])
μ_gompertz_treated_error = μ_gompertz_treated .± z .* sqrt.([FIM_G_treated[1,1], FIM_G_treated[2,2]])

# Distributions to be used for sample-drawing for Monte Carlo simulations
d_AI_control = Normal(μ_ageind_control[1], sqrt(FIM_AI_control[1]))
d_AI_treated = Normal(μ_ageind_treated[1], sqrt(FIM_AI_treated[1]))

d_G_control = MvNormal(μ_gompertz_control, FIM_G_control)
d_G_treated = MvNormal(μ_gompertz_treated, FIM_G_treated)

d_L_control = MvNormal(μ_logistic_control, convert(Matrix{Float64}, Hermitian(FIM_L_control)))
d_L_treated = MvNormal(μ_logistic_treated, convert(Matrix{Float64}, Hermitian(FIM_L_treated)))

# Transforming back for Logistic case
μ_logistic_control = exp.(μ_logistic_control)
μ_logistic_control_error = exp.(μ_logistic_control_error)
μ_logistic_treated = exp.(μ_logistic_treated)
μ_logistic_treated_error = exp.(μ_logistic_treated_error)

end
