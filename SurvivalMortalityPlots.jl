include("ParameterFitting.jl");

using .ParameterFitting # exported functions and parameters from ParameterFitting module. Includes: ages_control, ages_treated, Mortality, Survival1,  μ_ageind_control, μ_ageind_control_error, μ_ageind_treated, μ_ageind_treated_error, μ_logistic_control, μ_logistic_control_error, μ_logistic_treated, μ_logistic_treated_error, μ_gompertz_control, μ_gompertz_control_error, μ_gompertz_treated, μ_gompertz_treated_error, KM_control, KM_treated, CIU_C, CIU_T, CIL_C, CIL_T, d_AI_control, d_AI_treated, d_G_control, d_G_treated, d_L_control, d_L_treated

using Plots, LaTeXStrings, Random, Statistics

# Function to be used in sampling to get the same random sample
function rr(i)
    rng = MersenneTwister(i);
    return rng
end

ns = 10000 # number of simulations
CL = 0.95 # confidence for uncertainty from simulations
q = (1-CL)/2 # used for lower and upper (1-q) confidence level

μs_AI_control = rand(rr(123), d_AI_control, ns);
μs_AI_treated = rand(rr(123), d_AI_treated, ns);

μs_L_control = exp.(rand(rr(123), d_L_control, ns)); # transforming back
μs_L_treated = exp.(rand(rr(123), d_L_treated, ns)); # transforming back

μs_G_control = rand(rr(123), d_G_control, ns);
μs_G_treated = rand(rr(123), d_G_treated, ns);

# Useful parameters
a = collect(0:65); # array for ages 0-65

# Survival and Mortality functions
# Age-Independent
survival_AI_control = zeros(length(a), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    survival_AI_control[:, i] = Survival1(a, μ, "Age-Independent")
end

# lower and upper bounds
survival_AI_control_l = [quantile(survival_AI_control[i, :], q) for i in 1:length(a)]
survival_AI_control_u = [quantile(survival_AI_control[i, :], 1-q) for i in 1:length(a)]

survival_AI_treated = zeros(length(a), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    survival_AI_treated[:, i] = Survival1(a, μ, "Age-Independent")
end

survival_AI_treated_l = [quantile(survival_AI_treated[i, :], q) for i in 1:length(a)]
survival_AI_treated_u = [quantile(survival_AI_treated[i, :], 1-q) for i in 1:length(a)]

survival_control_ageind = Survival1(a, μ_ageind_control, "Age-Independent");
survival_treated_ageind = Survival1(a, μ_ageind_treated, "Age-Independent");

mortality_AI_control = zeros(length(a), ns)
for i in 1:ns
    μ = [μs_AI_control[i]]
    mortality_AI_control[:, i] = Mortality(a, μ, "Age-Independent")
end

# lower and upper bounds
mortality_AI_control_l = [quantile(mortality_AI_control[i, :], q) for i in 1:length(a)]
mortality_AI_control_u = [quantile(mortality_AI_control[i, :], 1-q) for i in 1:length(a)]

mortality_AI_treated = zeros(length(a), ns)
for i in 1:ns
    μ = [μs_AI_treated[i]]
    mortality_AI_treated[:, i] = Mortality(a, μ, "Age-Independent")
end

mortality_AI_treated_l = [quantile(mortality_AI_treated[i, :], q) for i in 1:length(a)]
mortality_AI_treated_u = [quantile(mortality_AI_treated[i, :], 1-q) for i in 1:length(a)]

mortality_control_ageind = Mortality(a, μ_ageind_control, "Age-Independent");
mortality_treated_ageind = Mortality(a, μ_ageind_treated, "Age-Independent");

# Survival fits with Kaplan-Meier estimator representing the data
survival_control_ageind_plot = plot(KM_control.times, CIL_C, fillrange=CIU_C, linetype=:steppost, ylims=(0,1), xlims=(0, 65), fillalpha=0.1, c=:orchid, l=0, show=true, minorgrid=false, label=false, tick_direction=:out);
plot!(KM_control.times, KM_control.survival, linetype=:steppost, label="Kaplan-Meier Estimator", c=:violet, l=2);
plot!(a, survival_AI_control_l, fillrange=survival_AI_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, survival_control_ageind, color=:dodgerblue, line=2, fillalpha=0.1, label="survival fit");
display(survival_control_ageind_plot)

survival_treated_ageind_plot = plot(KM_treated.times, CIL_T, fillrange=CIU_T, linetype=:steppost, ylims=(0,1), xlims=(0, 65), fillalpha=0.1, c=:orchid, l=0, show=true, minorgrid=false, label=false, tick_direction=:out);
plot!(KM_treated.times, KM_treated.survival, linetype=:steppost, label="", c=:violet, l=2);
plot!(a, survival_AI_treated_l, fillrange=survival_AI_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a, survival_treated_ageind, color=:chocolate4, line=2, fillalpha=0.1, label="");
display(survival_treated_ageind_plot)

# Mortality comparison
mortality_comparison_ageind = plot(a, mortality_AI_control_l, fillrange=mortality_AI_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, mortality_control_ageind, label="control", line=2, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Age-Ind)", show=true, ylims=(0.0,0.5), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
plot!(a, mortality_AI_treated_l, fillrange=mortality_AI_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a[1:33], mortality_treated_ageind[1:33], label="treated", line=2, color=:chocolate4);
plot!(a[33:end], mortality_treated_ageind[33:end], label="", line=(2, :dash), color=:chocolate4);
display(mortality_comparison_ageind)

# Logistic
survival_L_control = zeros(length(a), ns)
for i in 1:ns
    μ = μs_L_control[:, i]
    survival_L_control[:, i] = Survival1(a, μ, "Logistic")
end

# lower and upper bounds
survival_L_control_l = [quantile(survival_L_control[i, :], q) for i in 1:length(a)]
survival_L_control_u = [quantile(survival_L_control[i, :], 1-q) for i in 1:length(a)]

survival_L_treated = zeros(length(a), ns)
for i in 1:ns
    μ = μs_L_treated[:, i]
    survival_L_treated[:, i] = Survival1(a, μ, "Logistic")
end

survival_L_treated_l = [quantile(survival_L_treated[i, :], q) for i in 1:length(a)]
survival_L_treated_u = [quantile(survival_L_treated[i, :], 1-q) for i in 1:length(a)]

survival_control_logistic = Survival1(a, μ_logistic_control, "Logistic");
survival_treated_logistic = Survival1(a, μ_logistic_treated, "Logistic");

mortality_L_control = zeros(length(a), ns)
for i in 1:ns
    μ = μs_L_control[:, i]
    mortality_L_control[:, i] = Mortality(a, μ, "Logistic")
end

# lower and upper bounds
mortality_L_control_l = [quantile(mortality_L_control[i, :], q) for i in 1:length(a)]
mortality_L_control_u = [quantile(mortality_L_control[i, :], 1-q) for i in 1:length(a)]

mortality_L_treated = zeros(length(a), ns)
for i in 1:ns
    μ = μs_L_treated[:, i]
    mortality_L_treated[:, i] = Mortality(a, μ, "Logistic")
end

mortality_L_treated_l = [quantile(mortality_L_treated[i, :], q) for i in 1:length(a)]
mortality_L_treated_u = [quantile(mortality_L_treated[i, :], 1-q) for i in 1:length(a)]

mortality_control_logistic = Mortality(a, μ_logistic_control, "Logistic");
mortality_treated_logistic = Mortality(a, μ_logistic_treated, "Logistic");

# Survival fits
survival_control_logistic_plot = plot(KM_control.times, CIL_C, fillrange=CIU_C, linetype=:steppost, ylims=(0, 1), xlims =(0, 65), fillalpha = 0.1, c=:orchid, l=0, show=true, minorgrid=false, label=:false, tick_direction=:out);
plot!(KM_control.times, KM_control.survival, linetype=:steppost, label="Kaplan-Meier Estimator", c=:violet, l=2);
plot!(a, survival_L_control_l, fillrange=survival_L_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, survival_control_logistic, color=:dodgerblue, line=2, label="survival fit");
display(survival_control_logistic_plot)

survival_treated_logistic_plot = plot(KM_treated.times, CIL_T, fillrange=CIU_T, linetype=:steppost, ylims=(0, 1), xlims =(0, 65), fillalpha = 0.1, c=:orchid, l=0, show=true, minorgrid=false, label=:false, tick_direction=:out);
plot!(KM_treated.times, KM_treated.survival, linetype=:steppost, label="", c=:violet, l=2);
plot!(a, survival_L_treated_l, fillrange=survival_L_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a, survival_treated_logistic, color=:chocolate4, line=2, label="");
display(survival_treated_logistic_plot)

# Mortality comparison
mortality_comparison_logistic = plot(a, mortality_L_treated_l, fillrange=mortality_L_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a[1:33], mortality_treated_logistic[1:33], label="treated", line=2, color=:chocolate4);
plot!(a[33:end], mortality_treated_logistic[33:end], label="", line=(2, :dash), color=:chocolate4);
plot!(a, mortality_L_control_l, fillrange=mortality_L_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, mortality_control_logistic, label="control", line=2, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Logistic)", show=true, ylims=(0.0,0.65), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
display(mortality_comparison_logistic)

# Gompertz
survival_G_control = zeros(length(a), ns)
for i in 1:ns
    μ = μs_G_control[:, i]
    survival_G_control[:, i] = Survival1(a, μ, "Gompertz")
end

# lower and upper bounds
survival_G_control_l = [quantile(survival_G_control[i, :], q) for i in 1:length(a)]
survival_G_control_u = [quantile(survival_G_control[i, :], 1-q) for i in 1:length(a)]

survival_G_treated = zeros(length(a), ns)
for i in 1:ns
    μ = μs_G_treated[:, i]
    survival_G_treated[:, i] = Survival1(a, μ, "Gompertz")
end

survival_G_treated_l = [quantile(survival_G_treated[i, :], q) for i in 1:length(a)]
survival_G_treated_u = [quantile(survival_G_treated[i, :], 1-q) for i in 1:length(a)]

survival_control_gompertz = Survival1(a, μ_gompertz_control, "Gompertz");
survival_treated_gompertz = Survival1(a, μ_gompertz_treated, "Gompertz");

mortality_G_control = zeros(length(a), ns)
for i in 1:ns
    μ = μs_G_control[:, i]
    mortality_G_control[:, i] = Mortality(a, μ, "Gompertz")
end

# lower and upper bounds
mortality_G_control_l = [quantile(mortality_G_control[i, :], q) for i in 1:length(a)]
mortality_G_control_u = [quantile(mortality_G_control[i, :], 1-q) for i in 1:length(a)]

mortality_G_treated = zeros(length(a), ns)
for i in 1:ns
    μ = μs_G_treated[:, i]
    mortality_G_treated[:, i] = Mortality(a, μ, "Gompertz")
end

mortality_G_treated_l = [quantile(mortality_G_treated[i, :], q) for i in 1:length(a)]
mortality_G_treated_u = [quantile(mortality_G_treated[i, :], 1-q) for i in 1:length(a)]

mortality_control_gompertz = Mortality(a, μ_gompertz_control, "Gompertz");
mortality_treated_gompertz = Mortality(a, μ_gompertz_treated, "Gompertz");

# Survival fits with data
survival_control_gompertz_plot = plot(KM_control.times, CIL_C, fillrange=CIU_C, linetype=:steppost, ylims=(0, 1), xlims =(0, 65), fillalpha = 0.1, c=:orchid, l=0, show=true, minorgrid=false, label=:false, tick_direction=:out);
plot!(KM_control.times, KM_control.survival, linetype=:steppost, label="Kaplan-Meier Estimator", c=:violet, l=2);
plot!(a, survival_G_control_l, fillrange=survival_G_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, survival_control_gompertz, color=:dodgerblue, line=2, label="survival fit");
display(survival_control_gompertz_plot)

survival_treated_gompertz_plot = plot(KM_treated.times, CIL_T, fillrange=CIU_T, linetype=:steppost, ylims=(0, 1), xlims =(0, 65), fillalpha = 0.1, c=:orchid, l=0, show=true, minorgrid=false, label=:false, tick_direction=:out);
plot!(KM_treated.times, KM_treated.survival, linetype=:steppost, label="", c=:violet, l=2);
plot!(a, survival_G_treated_l, fillrange=survival_G_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a, survival_treated_gompertz, color=:chocolate4, line=2, label="");
display(survival_treated_gompertz_plot)


# Mortality comparison
mortality_comparison_gompertz = plot(a, mortality_G_treated_l, fillrange=mortality_G_treated_u, l=0, label=:false, fillalpha=0.1, c=:chocolate4);
plot!(a[1:33], mortality_treated_gompertz[1:33], label="treated", line=2, color=:chocolate4);
plot!(a[33:end], mortality_treated_gompertz[33:end], label="", line=(2, :dash), color=:chocolate4);
plot!(a, mortality_G_control_l, fillrange=mortality_G_control_u, l=0, label=:false, fillalpha=0.1, c=:dodgerblue);
plot!(a, mortality_control_gompertz, label="control", line=2, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Gompertz)", show=true, ylims=(0.0,0.5), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
display(mortality_comparison_gompertz)


# Figures for paper - plot panels
# Survival vs Control
# Age-Independent, Logistic, Gompertz
ageindcontrol = plot(survival_control_ageind_plot, legend=false);
ageindtreated = plot(survival_treated_ageind_plot, legendtitle="Age-Independent", foreground_color_legend=:false, legendtitlefontsize=20);
logisticcontrol = plot(survival_control_logistic_plot, legend=false, ylabel="Proportion of surviving mosquitoes");
logistictreated = plot(survival_treated_logistic_plot, legendtitle="Logistic", foreground_color_legend=:false, legendtitlefontsize=20);
gompertzcontrol = plot(survival_control_gompertz_plot, legend=false, xlabel="age (in days)");
gompertztreated = plot(survival_treated_gompertz_plot, legendtitle="Gompertz", foreground_color_legend=:false, legendtitlefontsize=20, xlabel="age (in days)");

survivalfits = [ageindcontrol, ageindtreated,
                logisticcontrol, logistictreated,
                gompertzcontrol, gompertztreated]

legendplotsurv = plot([], showaxis = false, grid = false, color=:orchid, label = "Kaplan-Meier Estimator", foreground_color_legend=nothing, legend=:outertop, l=2, ticks=false);
plot!([], showaxis = false, grid = false, color=:dodgerblue, line=2, label="control fits", legendfontsize=20, ticks=false);
plot!([], showaxis = false, grid = false, color=:chocolate4, line=2, label="treated fits", legendfontsize=20, ticks=false);

l = @layout [
    grid(3,2)
    a{0.05h}
];

SurvivalPlotsPanel = plot(survivalfits..., legendplotsurv, layout=l, tick_direction=:out, minorgrid=:false, size=(1750,1750), yguidefont=30, xtickfont=20, ytickfont=20, xguidefont=20, tickfontfamily="Times", font="Times", title=["Control" "Treated" "" "" "" "" ""], titlefontsize=30, show=true);
display(SurvivalPlotsPanel)

# Mortality comparisons
ageindmortality = plot(a, mortality_AI_control_l, fillrange=mortality_AI_control_u, fillalpha=0.1, label=:false, c=:dodgerblue, l=0);
plot!(a, mortality_control_ageind, label="", color=:dodgerblue, line=3, fillalpha=0.1);
plot!(a, mortality_AI_treated_l, fillrange=mortality_AI_treated_u, fillalpha=0.1, label=:false, c=:chocolate4, l=0);
plot!(a[1:33], mortality_treated_ageind[1:33], color=:chocolate4, line=3, label="", legendtitle="Age-Independent", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing);
plot!(a[33:end], mortality_treated_ageind[33:end], label="", line=(3, :dash), color=:chocolate4);

logisticmortality = plot(a, mortality_L_treated_l, fillrange=mortality_L_treated_u, fillalpha=0.1, label=:false, c=:chocolate4, l=0);
plot!(a[1:33], mortality_treated_logistic[1:33], color=:chocolate4, line=3, label="", legendtitle="Logistic", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing);
plot!(a[33:end], mortality_treated_logistic[33:end], label="", line=(3, :dash), color=:chocolate4, ylabel=L"\text{mortality rate (days^{-1})");
plot!(a, mortality_L_control_l, fillrange=mortality_L_control_u, fillalpha=0.1, label=:false, c=:dodgerblue, l=0);
plot!(a, mortality_control_logistic, label="", color=:dodgerblue, line=3, fillalpha=0.1);


gompertzmortality = plot(a, mortality_G_control_l, fillrange=mortality_G_control_u, fillalpha=0.1, label=:false, c=:dodgerblue, l=0);
plot!(a, mortality_control_gompertz, label="", color=:dodgerblue, line=3, fillalpha=0.1);
plot!(a, mortality_G_treated_l, fillrange=mortality_G_treated_u, fillalpha=0.1, label=:false, c=:chocolate4, l=0);
plot!(a, mortality_treated_gompertz, color=:chocolate4, line=3, label="", legendtitle="Gompertz", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing, xlabel="age (in days)");

legendplotmort = plot([], showaxis = false, grid = false, color=:dodgerblue, line=3, label = "control", legend=:outertop);
plot!([], showaxis = false, grid = false, color=:chocolate4, line=3, label="treated", legendfontsize=20);

mortalitycomparison = [ageindmortality, logisticmortality, gompertzmortality]

l = @layout [
    grid(3,1)
    m{0.05h}
]

MortalityPlotsPanel = plot(mortalitycomparison..., legendplotmort, layout=l, tick_direction=:out, minorgrid=:false, size=(1250,1750), xguidefont=20, yguidefont=30, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", ylims=(0, 0.6), xlims=(0, 64), xticks=(0:4:64), yticks=(0:0.1:0.6), show=true);
display(MortalityPlotsPanel)
