module ParameterFitting

export ages_control, ages_treated, proportion_alive_control, proportion_alive_treated, Mortality, Survival, μ_ageind_control, μ_ageind_control_error, μ_ageind_treated, μ_ageind_treated_error, μ_logistic_control, μ_logistic_control_error, μ_logistic_treated, μ_logistic_treated_error, μ_gompertz_control, μ_gompertz_control_error, μ_gompertz_treated, μ_gompertz_treated_error, fit_control_ageind, fit_treated_ageind, fit_control_logistic, fit_treated_logistic, fit_control_gompertz, fit_treated_gompertz

using XLSX, DataFrames, LsqFit, Measurements, Statistics

# Edit file paths accordingly

# Load survival data - CONTROL
df_survival_control = DataFrame(XLSX.readtable("../Supporting_information/S5_Table.xlsx", "SurvivalControl")...);

# Load survival data - TREATED
df_survival_treated = DataFrame(XLSX.readtable("../Supporting_information/S6_Table.xlsx", "SurvivalTreated")...);

# Ages and proportion alive from data including additional data-point
ages_control = vcat(0, df_survival_control.Age)
ages_treated = vcat(0, df_survival_treated.Age)
# Best-case scenario
proportion_alive_control = vcat(1.0, df_survival_control.Total/200)
proportion_alive_treated = vcat(1.0, df_survival_treated.Total/200)

# Uncomment for worst-case scenario
# proportion_alive_control = vcat(1.01, df_survival_control.Total/200)
# proportion_alive_treated = vcat(1.01, df_survival_treated.Total/200)

# Initial parameter values
μ₀_ageind = [0.0];
μ₀_logistic = [1.0, 0.1, 10.0];
μ₀_gompertz = [0.1, 0.1];

# Mortality function
function Mortality(a, μ, case)
    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            return μ .* ones(length(a))
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            return μ[1] ./ (1 .+ exp.(μ[2] .* (-a .+ μ[3])))
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            return μ[1] .* exp.(μ[2] .* a)
        end
    else
        return "Error: Invalid case"
    end
end

# Survival function
function Survival(a, μ, case)
    s0 = 1.00 # use this for best-case scenario
    # s0 = 1.01 # use this for worst-case scenario
    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            return  s0 * exp.(-μ.*a)
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            return s0 * (exp.(-((μ[1] .* (a .* μ[2] .- log.(1 .+ exp.(μ[2] * μ[3]))) .+ log.(1 .+ exp.(μ[2] .* (-a .+ μ[3]))))) ./ μ[2]))
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            return s0 * exp.(-μ[1] .* (exp.(μ[2] .* a) .- 1) ./ μ[2])
        end
    else
        return "Error: Invalid case"
    end
end

# Fitting
function fit(case, ages, proportion_alive, params_iv, lb, ub)
    if case != "Age-Independent" && case != "Logistic" && case != "Gompertz"
        return "Error: Invalid case"
    end
    if length(ages) != length(proportion_alive)
        return "Error: No data agreement (control/treated)"
    end
    if case == "Age-Independent" && length(params_iv) != 1
        return "Error: Invalid length of params_bounds"
    end
    if case == "Gompertz" && length(params_iv) != 2
        return "Error: Invalid length of params_bounds"
    end
    if case == "Logistic" && length(params_iv) != 3
        return "Error: Invalid length of params_bounds"
    end
    function model(a, μ)
        return Survival(a, μ, case)
    end
    fit = curve_fit(model, ages, proportion_alive, params_iv, lower=lb, upper=ub)
    return fit
end

# Parameter estimation from fit - breaking down the various elements obtained from the fit function
struct ParametersEstimation{μ, μ_error}
    param::μ
    error::μ_error
end

function parameter_estimation(fit)
    μ_params = coef(fit)
    μ_params_error = coef(fit) .± margin_error(fit)
    return ParametersEstimation(μ_params, μ_params_error)
end

# Fitting results
# Age-Independent
fit_control_ageind = fit("Age-Independent", ages_control, proportion_alive_control, μ₀_ageind, [0.0], [+Inf]);
fit_treated_ageind = fit("Age-Independent", ages_treated, proportion_alive_treated, μ₀_ageind, [0.0], [+Inf]);

μ_ageind_control = parameter_estimation(fit_control_ageind).param;
μ_ageind_control_error = parameter_estimation(fit_control_ageind).error;
μ_ageind_treated = parameter_estimation(fit_treated_ageind).param;
μ_ageind_treated_error = parameter_estimation(fit_treated_ageind).error;

# Logistic
fit_control_logistic = fit("Logistic", ages_control, proportion_alive_control, μ₀_logistic, [0.0, 0.0, 0.0], [+Inf, +Inf, +Inf]);
fit_treated_logistic = fit("Logistic", ages_treated, proportion_alive_treated, μ₀_logistic, [0.0, 0.0, 0.0], [5.0, 1.0, 40.0]);

μ_logistic_control = parameter_estimation(fit_control_logistic).param;
μ_logistic_control_error = parameter_estimation(fit_control_logistic).error;
μ_logistic_treated = parameter_estimation(fit_treated_logistic).param;
μ_logistic_treated_error = parameter_estimation(fit_treated_logistic).error;

# Gompertz
fit_control_gompertz = fit("Gompertz", ages_control, proportion_alive_control, μ₀_gompertz, [0.0, 0.0], [+Inf, +Inf]);
fit_treated_gompertz = fit("Gompertz", ages_treated, proportion_alive_treated, μ₀_gompertz, [0.0, 0.0], [+Inf, +Inf]);

μ_gompertz_control = parameter_estimation(fit_control_gompertz).param;
μ_gompertz_control_error = parameter_estimation(fit_control_gompertz).error;
μ_gompertz_treated = parameter_estimation(fit_treated_gompertz).param;
μ_gompertz_treated_error = parameter_estimation(fit_treated_gompertz).error;

end
