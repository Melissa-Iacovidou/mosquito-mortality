include("ParameterFitting.jl");

using .ParameterFitting # exported functions and parameters from ParameterFitting module. Includes: ages_control, ages_treated, proportion_alive_control, proportion_alive_treated, Mortality, Survival, μ_ageind_control, μ_ageind_control_error, μ_ageind_treated, μ_ageind_treated_error, μ_logistic_control, μ_logistic_control_error, μ_logistic_treated, μ_logistic_treated_error, μ_gompertz_control, μ_gompertz_control_error, μ_gompertz_treated, μ_gompertz_treated_error, fit_control_ageind, fit_treated_ageind, fit_control_logistic, fit_treated_logistic, fit_control_gompertz, fit_treated_gompertz

using Plots, LaTeXStrings, Measurements, Statistics

# Useful parameters
a = collect(0:65); # array for ages 0-65

# Standard Error from data
se_survival_control = std(proportion_alive_control) ./ sqrt.(proportion_alive_control * 200);
se_survival_treated = std(proportion_alive_treated) ./ sqrt.(proportion_alive_treated * 200);


# Survival and Mortality functions (with and without error)

# Age-Independent
# Survival - Contol
survival_control_ageind = Survival(a, μ_ageind_control, "Age-Independent");
survival_control_ageind_error = Survival(a, μ_ageind_control_error, "Age-Independent");
survival_control_ageind_unc = Measurements.uncertainty.(survival_control_ageind_error);

# Survival - Treated
survival_treated_ageind = Survival(a, μ_ageind_treated, "Age-Independent");
survival_treated_ageind_error = Survival(a, μ_ageind_treated_error, "Age-Independent");
survival_treated_ageind_unc = Measurements.uncertainty.(survival_treated_ageind_error);

# Mortality - Control
mortality_control_ageind = Mortality(a, μ_ageind_control, "Age-Independent");
mortality_control_ageind_error = Mortality(a, μ_ageind_control_error, "Age-Independent");
mortality_control_ageind_unc = Measurements.uncertainty.(mortality_control_ageind_error);

# Mortality - Treated
mortality_treated_ageind = Mortality(a, μ_ageind_treated, "Age-Independent");
mortality_treated_ageind_error = Mortality(a, μ_ageind_treated_error, "Age-Independent");
mortality_treated_ageind_unc = Measurements.uncertainty.(mortality_treated_ageind_error);

# Survival fits with data
survival_control_ageind_plot = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Age-Independent - Control", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(0:10:64));
plot!(a, survival_control_ageind, color=:dodgerblue, line=2, ribbon=survival_control_ageind_unc, fillalpha=0.1, label="survival fit");
display(survival_control_ageind_plot)

survival_treated_ageind_plot = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Age-Independent - Treated", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(4:10:64));
plot!(a, survival_treated_ageind, color=:chocolate4, line=2, ribbon=survival_treated_ageind_unc, fillalpha=0.1, label="survival fit");
display(survival_treated_ageind_plot)

# Mortality comparison
mortality_comparison_ageind = plot(a, mortality_control_ageind, label="control", line=3, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Age-Ind)", show=true, ribbon=mortality_control_ageind_unc, fillalpha=0.1, ylims=(0.0,0.5), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
plot!(a[1:33], mortality_treated_ageind[1:33], ribbon=mortality_treated_ageind_unc[1:33], fillalpha=0.1, label="treated", line=3, color=:chocolate4);
plot!(a[33:end], mortality_treated_ageind[33:end], ribbon=mortality_treated_ageind_unc[33:end], fillalpha=0.1, label="", line=(3, :dash), color=:chocolate4);
display(mortality_comparison_ageind)

# Logistic
# Survival - Control
survival_control_logistic = Survival(a, μ_logistic_control, "Logistic");
survival_control_logistic_error = Survival(a, μ_logistic_control_error, "Logistic");
survival_control_logistic_unc = Measurements.uncertainty.(survival_control_logistic_error);

# Survival - Treated
survival_treated_logistic = Survival(a, μ_logistic_treated, "Logistic");
survival_treated_logistic_error = Survival(a, μ_logistic_treated_error, "Logistic");
survival_treated_logistic_unc = Measurements.uncertainty.(survival_treated_logistic_error);

# Mortality - Control
mortality_control_logistic = Mortality(a, μ_logistic_control, "Logistic");
mortality_control_logistic_error = Mortality(a, μ_logistic_control_error, "Logistic");
mortality_control_logistic_unc = Measurements.uncertainty.(mortality_control_logistic_error);

# Mortality - Treated
mortality_treated_logistic = Mortality(a, μ_logistic_treated, "Logistic");
mortality_treated_logistic_error = Mortality(a, μ_logistic_treated_error, "Logistic");
mortality_treated_logistic_unc = Measurements.uncertainty.(mortality_treated_logistic_error);

# Survival fits with data
survival_control_logistic_plot = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Logistic - Control", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(0:4:64));
plot!(a, survival_control_logistic, color=:dodgerblue, line=2, ribbon=survival_control_logistic_unc, fillalpha=0.1, label="survival fit");
display(survival_control_logistic_plot)

survival_treated_logistic_plot = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Logistic - Treated", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(0:10:64));
plot!(a, survival_treated_logistic, color=:chocolate4, line=2, ribbon=survival_treated_logistic_unc, fillalpha=0.1, label="survival fit");
display(survival_treated_logistic_plot)


# Mortality comparison
mortality_comparison_logistic = plot(a, mortality_control_logistic, label="control", line=3, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Logistic)", show=true, ribbon=mortality_control_logistic_unc, fillalpha=0.1, ylims=(0.0,1.0), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
plot!(a[1:33], mortality_treated_logistic[1:33], ribbon=mortality_treated_logistic_unc[1:34], fillalpha=0.1, label="treated", line=3, color=:chocolate4);
plot!(a[33:end], mortality_treated_logistic[33:end], ribbon=mortality_treated_logistic_unc[33:end], fillalpha=0.1, label="", line=(3, :dash), color=:chocolate4);
display(mortality_comparison_logistic)

# Gompertz
# Survival - Control
survival_control_gompertz = Survival(a, μ_gompertz_control, "Gompertz");
survival_control_gompertz_error = Survival(a, μ_gompertz_control_error, "Gompertz");
survival_control_gompertz_unc = Measurements.uncertainty.(survival_control_gompertz_error);

# Survival - Treated
survival_treated_gompertz = Survival(a, μ_gompertz_treated, "Gompertz");
survival_treated_gompertz_error = Survival(a, μ_gompertz_treated_error, "Gompertz");
survival_treated_gompertz_unc = Measurements.uncertainty.(survival_treated_gompertz_error);

# Mortality = Control
mortality_control_gompertz = Mortality(a, μ_gompertz_control, "Gompertz");
mortality_control_gompertz_error = Mortality(a, μ_gompertz_control_error, "Gompertz");
mortality_control_gompertz_unc = Measurements.uncertainty.(mortality_control_gompertz_error);

# Mortality - Treated
mortality_treated_gompertz = Mortality(a, μ_gompertz_treated, "Gompertz");
mortality_treated_gompertz_error = Mortality(a, μ_gompertz_treated_error, "Gompertz");
mortality_treated_gompertz_unc = Measurements.uncertainty.(mortality_treated_gompertz_error);

# Survival fits with data
survival_control_gompertz_plot = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Gompertz - Control", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(0:10:64));
plot!(a, survival_control_gompertz, color=:dodgerblue, line=2, ribbon=survival_control_gompertz_unc, fillalpha=0.1, label="survival fit");
display(survival_control_gompertz_plot)

survival_treated_gompertz_plot = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01, 1.05), xlim=(-0.5, 64.5), show=true, minorgrid=false, yaxis="Proportion of alive mosquitoes", title="Survival function - Gompertz - Treated", color=:violet, msc=:orchid, msw=:0.5, label="proportion alive (data)", tick_direction=:out, xaxis="age (in days)", fontfamily="Times", xticks=(0:10:64));
plot!(a, survival_treated_gompertz, color=:chocolate4, line=2, ribbon=survival_treated_gompertz_unc, fillalpha=0.1, label="survival fit");
display(survival_treated_gompertz_plot)


# Mortality comparison
mortality_comparison_gompertz = plot(a, mortality_control_gompertz, label="control", line=3, color=:dodgerblue, xaxis="age (in days)", fontfamily="Times", yaxis=L"\textrm{mortality rate (days}^{-1})", title="Comparison of mortality for the two treatments (Gompertz)", show=true, ribbon=mortality_control_gompertz_unc, fillalpha=0.1, ylims=(0.0,1), tickdirection=:out, xticks=(4:6:64), xlims=(4,64));
plot!(a[1:33], mortality_treated_gompertz[1:33], ribbon=mortality_treated_gompertz_unc[1:34], fillalpha=0.1, label="treated", line=3, color=:chocolate4);
plot!(a[33:end], mortality_treated_gompertz[33:end], ribbon=mortality_treated_gompertz_unc[33:end], fillalpha=0.1, label="", line=(3, :dash), color=:chocolate4);
display(mortality_comparison_gompertz)


# Figures for paper - plot panels
# Survival vs Control
# Age-Independent, Logistic, Gompertz
ageindcontrol = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false);
plot!(a, survival_control_ageind, color=:dodgerblue, line=2, ribbon=survival_control_ageind_unc, fillalpha=0.1, label=false);

ageindtreated = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false, ytickfontcolor=:white);
plot!(a, survival_treated_ageind, color=:chocolate4, line=2, label=:false, ribbon=survival_treated_ageind_unc, fillalpha=0.1, legendtitle="Age-Independent", legend=:topright, foreground_color_legend=:false, legendtitlefontsize=20);

logisticcontrol = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false);
plot!(a, survival_control_logistic, color=:dodgerblue, line=2, label=false, ribbon=survival_control_logistic_unc, fillalpha=0.1, yaxis="Proportion of surviving mosquitoes");

logistictreated = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false, ytickfontcolor=:white);
plot!(a, survival_treated_logistic, label=false, color=:chocolate4, line=2, ribbon=survival_treated_logistic_unc, fillalpha=0.1, legendtitle="Logistic", foreground_color_legend=nothing, legendtitlefontsize=20);

gompertzcontrol = scatter(ages_control, proportion_alive_control .± se_survival_control, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false);
plot!(a, survival_control_gompertz, color=:dodgerblue, line=2, label=false, ribbon=survival_control_gompertz_unc, fillalpha=0.1, xaxis="age (in days)", xguidefont=20);

gompertztreated = scatter(ages_treated, proportion_alive_treated .± se_survival_treated, ylim=(-0.01,1.05), xlim=(-0.5,64.5), color=:violet, msw=0.5, msc=:orchid, xticks=(0:10:64), label=:false, ytickfontcolor=:white);
plot!(a, survival_treated_gompertz, label=false, color=:chocolate4, line=2, ribbon=survival_treated_gompertz_unc, fillalpha=0.1, xaxis="age (in days)", legendtitle="Gompertz", foreground_color_legend=nothing, legendtitlefontsize=20, xguidefont=20);

survivalfits = [ageindcontrol, ageindtreated,
                logisticcontrol, logistictreated,
                gompertzcontrol, gompertztreated]

controltitle = plot([], xlabel="Control", grid=false, axis=false, leg=false, xguidefont=30);
treatedtitle = plot([], xlabel="Treated", grid=false, axis=false, leg=false, xguidefont=30);

legendplotsurv = scatter([], showaxis = false, grid = false, color=:violet, label = "proportion alive", foreground_color_legend=nothing, legend=:outertop, msw=0.5, msc=:orchid);
plot!([], showaxis = false, grid = false, color=:dodgerblue, line=2, label="control fits", legendfontsize=20);
plot!([], showaxis = false, grid = false, color=:chocolate4, line=2, label="treated fits", legendfontsize=20);

l = @layout [
    grid(1,2){0.001h}
    grid(3,2)
    a{0.05h}
];

SurvivalPlotsPanel = plot(controltitle, treatedtitle, survivalfits..., legendplotsurv, layout=l, tick_direction=:out, minorgrid=:false, size=(1750,1750), yguidefont=30, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true);
display(SurvivalPlotsPanel)


# Mortality comparisons
ageindmortality = plot(a, mortality_control_ageind, label="", ribbon=mortality_control_ageind_unc, color=:dodgerblue, line=3, fillalpha=0.1);
plot!(a[1:33], mortality_treated_ageind[1:33], ribbon=mortality_treated_ageind_unc, color=:chocolate4, line=3, label="", fillalpha=0.1, legendtitle="Age-Independent", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing);
plot!(a[33:end], mortality_treated_ageind[33:end], ribbon=mortality_treated_ageind_unc[33:end], fillalpha=0.1, label="", line=(3, :dash), color=:chocolate4);

logisticmortality = plot(a, mortality_control_logistic, line=3, color=:dodgerblue, yaxis=L"\textrm{mortality rate (days}^{-1})", ribbon=mortality_control_logistic_unc, fillalpha=0.1, label="");
plot!(a[1:33], mortality_treated_logistic[1:33], line=3, color=:chocolate4, label="", legendtitle="Logistic", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing, ribbon=mortality_treated_logistic_unc, fillalpha=0.1);
plot!(a[33:end], mortality_treated_logistic[33:end], label="", line=(:dash, 3), color=:chocolate4, ribbon=mortality_treated_logistic_unc[33:end], fillalpha=0.1);

gompertzmortality = plot(a, mortality_control_gompertz, line=3, color=:dodgerblue, ribbon=mortality_control_gompertz_unc, fillalpha=0.1, label="");
plot!(a, mortality_treated_gompertz, line=3, color=:chocolate4, label="", legendtitle="Logistic", legend=:topleft, legendtitlefontsize=20, foreground_color_legend=nothing, ribbon=mortality_treated_gompertz_unc, fillalpha=0.1);

titleplot = plot([], title="Comparison of mortality for the two treatments", grid=false, axis=false, leg=false, titlefont=30);

legendplotmort = plot([], showaxis = false, grid = false, color=:dodgerblue, line=3, label = "control", legend=:outertop);
plot!([], showaxis = false, grid = false, color=:chocolate4, line=3, label="treated", legendfontsize=20);

mortalitycomparison = [ageindmortality, logisticmortality, gompertzmortality]

l = @layout [
    a{0.001h}
    grid(3,1)
    m{0.05h}
]

MortalityPlotsPanel = plot(titleplot, mortalitycomparison..., legendplotmort, layout=l, tick_direction=:out, minorgrid=:false, size=(1250,1750), xguidefont=20, yguidefont=30, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", ylims=(0,1), xlims=(0, 64), xticks=(0:4:64), show=true);
display(MortalityPlotsPanel)
