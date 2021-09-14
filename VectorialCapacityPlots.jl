include("ParameterFitting.jl")
include("VCCalculations.jl")

using .ParameterFitting, .VCCalculations

using Plots, LaTeXStrings, Measurements, LsqFit

# Useful parameters
α = 0.25; # bite rate
σ = 0.09698706953718507; # incubation rate for EIP distribution - obtained from the ErlangEIP.nb mathematica file
λ = 31*σ; # for Erlang distribution with shape parameter k = 31
bites_k = collect(0:50); # bites in step 2
a₀_survEIP = collect(0:64); # a₀ in step 1
a₁_bites = collect(5:35); # a₁ in step 2
a₀_bites = collect(0:35); # a₀ in step 3


# Step 1: What is the probability of a mosquito surviving the EIP given it takes an infectious blood-meal at age a₀?

surviving_EIP_control_ageind = surviving_EIP(μ_ageind_control, a₀_survEIP, "Age-Independent", "Erlang", σ, λ);
surviving_EIP_control_ageind_error = surviving_EIP(μ_ageind_control_error, a₀_survEIP, "Age-Independent", "Erlang", σ, λ);
surviving_EIP_control_ageind_unc = Measurements.uncertainty.(surviving_EIP_control_ageind_error);

surviving_EIP_control_logistic = surviving_EIP(μ_logistic_control, a₀_survEIP, "Logistic", "Erlang", σ, λ);
surviving_EIP_control_logistic_error = surviving_EIP(μ_logistic_control_error, a₀_survEIP, "Logistic", "Erlang", σ, λ);
surviving_EIP_control_logistic_unc = Measurements.uncertainty.(surviving_EIP_control_logistic_error);

surviving_EIP_control_gompertz = surviving_EIP(μ_gompertz_control, a₀_survEIP, "Gompertz", "Erlang", σ, λ);
surviving_EIP_control_gompertz_error = surviving_EIP(μ_gompertz_control_error, a₀_survEIP, "Gompertz", "Erlang", σ, λ);
surviving_EIP_control_gompertz_unc = Measurements.uncertainty.(surviving_EIP_control_gompertz_error);

# Plot for control
surv_EIP_control = plot(a₀_survEIP, surviving_EIP_control_ageind, ribbon=surviving_EIP_control_ageind_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:purple4, tickdirection=:out, label="Age-Independent", fontfamily="Times", title="Surviving EIP control");
plot!(a₀_survEIP, surviving_EIP_control_logistic, ribbon=surviving_EIP_control_logistic_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:teal, tickdirection=:out, label="Logistic", fontfamily="Times");
plot!(a₀_survEIP, surviving_EIP_control_gompertz, ribbon=surviving_EIP_control_gompertz_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:orange, tickdirection=:out, label="Gompertz", fontfamily="Times");
display(surv_EIP_control)

surviving_EIP_treated_ageind = surviving_EIP(μ_ageind_treated, a₀_survEIP, "Age-Independent", "Erlang", σ, λ);
surviving_EIP_treated_ageind_error = surviving_EIP(μ_ageind_treated_error, a₀_survEIP, "Age-Independent", "Erlang", σ, λ);
surviving_EIP_treated_ageind_unc = Measurements.uncertainty.(surviving_EIP_treated_ageind_error);

surviving_EIP_treated_logistic = surviving_EIP(μ_logistic_treated, a₀_survEIP, "Logistic", "Erlang", σ, λ);
surviving_EIP_treated_logistic_error = surviving_EIP(μ_logistic_treated_error, a₀_survEIP, "Logistic", "Erlang", σ, λ);
surviving_EIP_treated_logistic_unc = Measurements.uncertainty.(surviving_EIP_treated_logistic_error);

surviving_EIP_treated_gompertz = surviving_EIP(μ_gompertz_treated, a₀_survEIP, "Gompertz", "Erlang", σ, λ);
surviving_EIP_treated_gompertz_error = surviving_EIP(μ_gompertz_treated_error, a₀_survEIP, "Gompertz", "Erlang", σ, λ);
surviving_EIP_treated_gompertz_unc = Measurements.uncertainty.(surviving_EIP_treated_gompertz_error);

# Plot for treated
surv_EIP_treated = plot(a₀_survEIP, surviving_EIP_treated_ageind, ribbon=surviving_EIP_treated_ageind_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:purple4, tickdirection=:out, label="Age-Independent", fontfamily="Times", title="Surviving EIP treated");
plot!(a₀_survEIP, surviving_EIP_treated_logistic, ribbon=surviving_EIP_treated_logistic_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:teal, tickdirection=:out, label="Logistic", fontfamily="Times");
plot!(a₀_survEIP, surviving_EIP_treated_gompertz, ribbon=surviving_EIP_treated_gompertz_unc, fillalpha=0.1, ylims=(0,1), xlim=(1,65), xticks=(1:8:65), line=2, c=:orange, tickdirection=:out, label="Gompertz", fontfamily="Times");
display(surv_EIP_treated)


# Step 2: How many bites will the mosquito take if it has survived the EIP?

# Control
bites_surv_EIP_control_ageind = no_bites(μ_ageind_control, a₁_bites, "Age-Independent").bites;
bites_surv_EIP_control_ageind_error = no_bites(μ_ageind_control_error, a₁_bites, "Age-Independent").bites;
means_bites_control_ageind = no_bites(μ_ageind_control_error, a₁_bites, "Age-Independent").means;

bites_surv_EIP_control_logistic = no_bites(μ_logistic_control, a₁_bites, "Logistic").bites;
bites_surv_EIP_control_logistic_error = no_bites(μ_logistic_control_error, a₁_bites, "Logistic").bites;
means_bites_control_logistic = no_bites(μ_logistic_control_error, a₁_bites, "Logistic").means;

bites_surv_EIP_control_gompertz = no_bites(μ_gompertz_control, a₁_bites, "Gompertz").bites;
bites_surv_EIP_control_gompertz_error = no_bites(μ_gompertz_control_error, a₁_bites, "Gompertz").bites;
means_bites_control_gompertz = no_bites(μ_gompertz_control_error, a₁_bites, "Gompertz").means;

heatmap(a₁_bites, bites_k, bites_surv_EIP_control_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Age-Independent Control");
scatter!(a₁_bites, means_bites_control_ageind, color=:white, marker=:square, label="", fillalpha=0.75)

heatmap(a₁_bites, bites_k, bites_surv_EIP_control_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Logistic Control");
scatter!(a₁_bites, means_bites_control_logistic, color=:white, marker=:square, label="", fillalpha=0.75)

heatmap(a₁_bites, bites_k, bites_surv_EIP_control_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Gompertz Control");
scatter!(a₁_bites, means_bites_control_gompertz, color=:white, marker=:square, label="", fillalpha=0.75)

# Treated
bites_surv_EIP_treated_ageind = no_bites(μ_ageind_treated, a₁_bites, "Age-Independent").bites;
bites_surv_EIP_treated_ageind_error = no_bites(μ_ageind_treated_error, a₁_bites, "Age-Independent").bites;
means_bites_treated_ageind = no_bites(μ_ageind_treated_error, a₁_bites, "Age-Independent").means;

bites_surv_EIP_treated_logistic = no_bites(μ_logistic_treated, a₁_bites, "Logistic").bites;
bites_surv_EIP_treated_logistic_error = no_bites(μ_logistic_treated_error, a₁_bites, "Logistic").bites;
means_bites_treated_logistic = no_bites(μ_logistic_treated_error, a₁_bites, "Logistic").means;

bites_surv_EIP_treated_gompertz = no_bites(μ_gompertz_treated, a₁_bites, "Gompertz").bites;
bites_surv_EIP_treated_gompertz_error = no_bites(μ_gompertz_treated_error, a₁_bites, "Gompertz").bites;
means_bites_treated_gompertz = no_bites(μ_gompertz_treated_error, a₁_bites, "Gompertz").means;

heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Age-Independent Treated");
scatter!(a₁_bites, means_bites_treated_ageind, color=:white, marker=:square, label="", fillalpha=0.75)

heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Logistic Treated");
scatter!(a₁_bites, means_bites_treated_logistic, color=:white, marker=:square, label="", fillalpha=0.75)

heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), yticks=0:3:12, fontfamily="Times", title="Gompertz Treated");
scatter!(a₁_bites, means_bites_treated_gompertz, color=:white, marker=:square, label="", fillalpha=0.75)


# Step 3: What is the expected number of bites a mosquito takes in its lifetime if it has taken an infectious blood-meal at age a₀?

# Control
bites_age_control_ageind = bites_age(μ_ageind_control, a₀_bites, "Age-Independent", "Erlang", σ, α, λ);
bites_age_control_ageind_error = bites_age(μ_ageind_control_error, a₀_bites, "Age-Independent", "Erlang", σ, α, λ);
bites_age_control_ageind_unc = Measurements.uncertainty.(bites_age_control_ageind_error);

bites_age_control_logistic = bites_age(μ_logistic_control, a₀_bites, "Logistic", "Erlang", σ, α, λ);
bites_age_control_logistic_error = bites_age(μ_logistic_control_error, a₀_bites, "Logistic", "Erlang", σ, α, λ);
bites_age_control_logistic_unc = Measurements.uncertainty.(bites_age_control_logistic_error);

bites_age_control_gompertz = bites_age(μ_gompertz_control, a₀_bites, "Gompertz", "Erlang", σ, α, λ);
bites_age_control_gompertz_error = bites_age(μ_gompertz_control_error, a₀_bites, "Gompertz", "Erlang", σ, α, λ);
bites_age_control_gompertz_unc = Measurements.uncertainty.(bites_age_control_gompertz_error);

bites_a0_control = plot(a₀_bites, bites_age_control_ageind, ribbon=bites_age_control_ageind_unc, fillalpha=0.1, ylims=(0,6.1), xlim=(0,35), xticks=(0:5:35), line=2, c=:purple4, tickdirection=:out, label="Age-Independent", fontfamily="Times", title=L"\textrm{Bites depending on}~a_0~\textrm{Control}", legend=:right);
plot!(a₀_bites, bites_age_control_logistic, ribbon=bites_age_control_logistic_unc, fillalpha=0.1, line=2, c=:teal, tickdirection=:out, label="Logistic", fontfamily="Times");
plot!(a₀_bites, bites_age_control_gompertz, ribbon=bites_age_control_gompertz_unc, fillalpha=0.1, line=2, c=:orange, tickdirection=:out, label="Gompertz", fontfamily="Times");
display(bites_a0_control)

# Treated
bites_age_treated_ageind = bites_age(μ_ageind_treated, a₀_bites, "Age-Independent", "Erlang", σ, α, λ);
bites_age_treated_ageind_error = bites_age(μ_ageind_treated_error, a₀_bites, "Age-Independent", "Erlang", σ, α, λ);
bites_age_treated_ageind_unc = Measurements.uncertainty.(bites_age_treated_ageind_error);

bites_age_treated_logistic = bites_age(μ_logistic_treated, a₀_bites, "Logistic", "Erlang", σ, α, λ);
bites_age_treated_logistic_error = bites_age(μ_logistic_treated_error, a₀_bites, "Logistic", "Erlang", σ, α, λ);
bites_age_treated_logistic_unc = Measurements.uncertainty.(bites_age_treated_logistic_error);

bites_age_treated_gompertz = bites_age(μ_gompertz_treated, a₀_bites, "Gompertz", "Erlang", σ, α, λ);
bites_age_treated_gompertz_error = bites_age(μ_gompertz_treated_error, a₀_bites, "Gompertz", "Erlang", σ, α, λ);
bites_age_treated_gompertz_unc = Measurements.uncertainty.(bites_age_treated_gompertz_error);

bites_a0_treated = plot(a₀_bites, bites_age_treated_ageind, ribbon=bites_age_treated_ageind_unc, fillalpha=0.1, ylims=(0,6.1), xlim=(0,35), xticks=(0:5:35), line=2, c=:purple4, tickdirection=:out, label="Age-Independent", fontfamily="Times", title=L"\textrm{Bites depending on}~a_0~\textrm{Treated}");
plot!(a₀_bites, bites_age_treated_logistic, ribbon=bites_age_treated_logistic_unc, fillalpha=0.1, line=2, c=:teal, tickdirection=:out, label="Logistic", fontfamily="Times");
plot!(a₀_bites, bites_age_treated_gompertz, ribbon=bites_age_treated_gompertz_unc, fillalpha=0.1, line=2, c=:orange, tickdirection=:out, label="Gompertz", fontfamily="Times");
display(bites_a0_treated)


# Step 4: What is the expected number of (infectious) bites a mosquito takes in its lifetime?

# Control
expected_bites_ageind_control = expected_bites(μ_ageind_control_error, "Age-Independent", "Erlang", σ, α, λ)

expected_bites_logistic_control = expected_bites(μ_logistic_control_error, "Logistic", "Erlang", σ, α, λ)

expected_bites_gompertz_control = expected_bites(μ_gompertz_control_error, "Gompertz", "Erlang", σ, α, λ)

# Treated
expected_bites_ageind_treated = expected_bites(μ_ageind_treated_error, "Age-Independent", "Erlang", σ, α, λ)

expected_bites_logistic_treated = expected_bites(μ_logistic_treated_error, "Logistic", "Erlang", σ, α, λ)

expected_bites_gompertz_treated = expected_bites(μ_gompertz_treated_error, "Gompertz", "Erlang", σ, α, λ)

# Vectorial Capacity (using the constants from VCCalculations.jl)
# Control
vc_control_ageind = VectorialCapacity(expected_bites_ageind_control)

vc_control_logistic = VectorialCapacity(expected_bites_logistic_control)

vc_control_gompertz = VectorialCapacity(expected_bites_gompertz_control)

# Treated
vc_treated_ageind = VectorialCapacity(expected_bites_ageind_treated)

vc_treated_logistic = VectorialCapacity(expected_bites_logistic_treated)

vc_treated_gompertz = VectorialCapacity(expected_bites_gompertz_treated)

# Decrease in VC between control and treatment
ageinderl = TreatmentDecrease(vc_control_ageind, vc_treated_ageind)

logisticerl = TreatmentDecrease(vc_control_logistic, vc_treated_logistic)

gompertzerl = TreatmentDecrease(vc_control_gompertz, vc_treated_gompertz)


# Plots for paper

probsurvEIPcontrol = plot(a₀_survEIP, surviving_EIP_control_ageind, legend=:topright, line=(3, :dot), xlabel=L"a_0~\textrm{(age infectious blood-meal taken)}", label=:false, color=:purple4, ylabel=L"P(\textrm{surviving EIP}~|~\textrm{blood-meal at}~a_0)", ribbon=surviving_EIP_control_ageind_unc, fillalpha=0.1, legendtitle="Control", foreground_color_legend=nothing, legendtitlefontsize=15, show=true);
plot!(a₀_survEIP, surviving_EIP_control_logistic, ribbon=surviving_EIP_control_logistic_unc, line=(3, :dash), label=:false, color=:teal, fillalpha=0.1);
plot!(a₀_survEIP, surviving_EIP_control_gompertz, ribbon=surviving_EIP_control_gompertz_unc, line=3, label=:false, color=:orange, fillalpha=0.1);

probsurvEIPtreated = plot(a₀_survEIP, surviving_EIP_treated_ageind, legend=:topright, line=(3, :dot), xlabel=L"a_0~\textrm{(age infectious blood-meal taken in days)}", label=:false, color=:purple4, ylabel="", ribbon=surviving_EIP_treated_ageind_unc, fillalpha=0.1, legendtitle="Treated", foreground_color_legend=nothing, legendtitlefontsize=15, show=true);
plot!(a₀_survEIP, surviving_EIP_treated_logistic, ribbon=surviving_EIP_treated_logistic_unc, line=(3, :dash), label=:false, color=:teal, fillalpha=0.1);
plot!(a₀_survEIP, surviving_EIP_treated_gompertz, ribbon=surviving_EIP_treated_gompertz_unc, line=3, label=:false, color=:orange, fillalpha=0.1);

probsurvEIP = [probsurvEIPcontrol, probsurvEIPtreated]

titleplot = plot([], grid=false, axis=false, leg=false, titlefont=20, title=L"\textrm{Probability of surviving EIP}~|~\textrm{infectious blood-meal at age} ~a_0");

legendplot = plot([], label="Gompertz", color=:orange, grid=false, showaxis=false, line=3, legend=:outertop, legendfontsize=10);
plot!([], label="Age-Independent", color=:purple4, line=(3, :dot));
plot!([], label="Logistic", color=:teal, line=(3, :dash), foreground_color_legend=nothing);

l = @layout[a{0.001h}
            grid(1, 2)
            b{0.1h}]
step1 = plot(titleplot, probsurvEIP..., legendplot, layout=l, ylims=(0,1.0), xlims=(0, 64), xticks=(0:8:64), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10);
display(step1)

# Step 2
mean_ageind_control = means_bites_control_ageind[11];
mean_ageind_treated = means_bites_treated_ageind[11];
mean_logistic_control = means_bites_control_logistic[11];
mean_logistic_treated = means_bites_treated_logistic[11];
mean_gompertz_control = means_bites_control_gompertz[11];
mean_gompertz_treated = means_bites_treated_gompertz[11];

pmf_ageind = scatter(bites_k, bites_surv_EIP_control_ageind_error[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), label="control - $mean_ageind_control");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_ageind_error[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, label="treated - $mean_ageind_treated", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=n~|\textrm{~exits EIP at~}a_1 = 15)");

pmf_logistic = scatter(bites_k, bites_surv_EIP_control_logistic_error[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), label="control - $mean_logistic_control");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_logistic_error[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, label="treated - $mean_logistic_treated", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=n~|\textrm{~exits EIP at~}a_1 = 15)");

pmf_gompertz = scatter(bites_k, bites_surv_EIP_control_gompertz_error[:,11], xlims=(-0.5, 20.5), ms=5, m=:rect, c=:dodgerblue, ylims=(0, 0.5), label="control - $mean_gompertz_control");
scatter!(collect(0.1:50.1), bites_surv_EIP_treated_gompertz_error[:,11], m=:rect, ms=5, tickdirection=:out, c=:chocolate4, label="treated - $mean_gompertz_treated", fontfamily="Times");
xticks!(0:1:50, ylabel=L"P(z=n~|\textrm{~exits EIP at~}a_1 = 15)", xlabel=L"\textrm{Number of bites}~(n)");

pmfs = [pmf_ageind, pmf_logistic, pmf_gompertz];

pmfs_panel = plot(pmfs..., layout=grid(3,1), fontfamily="Times", size=(600, 1000), show=true, tickfontsize=10, xguidefontsize=12.5);
display(pmfs_panel)


# Heatmap panel
heatmap_ageind_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), ylabel="Control", legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_ageind, color=:white, marker=:square, label="", alpha=0.75);

heatmap_ageind_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_ageind, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), ylabel="Treated", legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_ageind, color=:white, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}");

heatmap_logistic_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_logistic, color=:white, marker=:square, label="", alpha=0.75);

heatmap_logistic_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_logistic, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_logistic, color=:white, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}");

heatmap_gompertz_control = heatmap(a₁_bites, bites_k, bites_surv_EIP_control_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_control_gompertz, color=:white, marker=:square, label="", alpha=0.75);

heatmap_gompertz_treated = heatmap(a₁_bites, bites_k, bites_surv_EIP_treated_gompertz, clims=(0,0.5), xlim=(4.5, 35.5), ylim=(0,12), legend=:none, show=true, yguidefontsize=20);
scatter!(a₁_bites, means_bites_treated_gompertz, color=:white, marker=:square, label="", alpha=0.75, xlabel=L"a_1~\textrm{(age mosquito exits EIP in days)}");

heatmaps = [heatmap_ageind_control, heatmap_logistic_control, heatmap_gompertz_control, heatmap_ageind_treated, heatmap_logistic_treated, heatmap_gompertz_treated]

ageindtitle = plot([], title="Age-Independent", grid=false, axis=false, leg=false, titlefont=30);
logistictitle = plot([], title="Logistic", grid=false ,axis=false, leg=false, titlefont=30);
gompertztitle = plot([], title="Gompertz", grid=false ,axis=false, leg=false, titlefont=30);

yaxisplot = plot([], grid=false, axis=false, leg=false, ylabel="Number of bites (n)", mirror=:true, show=true, yguidefontsize=30);

colourbarplot = heatmap((0:0.001:0.5).*ones(501,1), legend=:none, xticks=:none, yticks=(1:100:501, string.(0:0.1:0.5)), ylabel=L"P(z=n~|~\textrm{exits EIP at}~a_1)", mirror=:true, yguidefontsize=30);

emptyspace = plot([], grid=false, axis=false, leg=false);

legendplot = scatter([], showaxis = false, grid = false, color=:white, alpha=0.75, marker=:square, label="mean no. of bites", foreground_color_legend=nothing, legend=:outertop, legendfontsize=10);

l = @layout[grid(1, 3){0.001h}
            a{0.001w} b{0.01w} grid(2,3) c{0.05w} d{0.001w}
            e{0.05h}]

heatmap_panel = plot(ageindtitle, logistictitle, gompertztitle, yaxisplot, emptyspace, heatmaps..., colourbarplot, emptyspace, legendplot, layout=l, fontfamily="Times", size=(1900, 1300), show=true, tickfontsize=15, xguidefontsize=20, yticks=0:3:12);
display(heatmap_panel)


# Step 3 plots

expected_bites_control = plot(a₀_bites, bites_age_control_ageind, ribbon=bites_age_control_ageind_unc, fillalpha=0.1, line=(3, :dot), color=:purple4, legendtitle="Control", legend=:topright, label="", foreground_color_legend=nothing, xlabel=L"a_0~\textrm{(age infectious blood-meal taken in days)}", ylabel=L"E(z~|~\textrm{inf. blood-meal at}~a_0)");
plot!(a₀_bites, bites_age_control_logistic, ribbon=bites_age_control_logistic_unc, fillalpha=0.1, line=(3, :dash), color=:teal, label="", foreground_color_legend=nothing);
plot!(a₀_bites, bites_age_control_gompertz, ribbon=bites_age_control_gompertz_unc, fillalpha=0.1, line=3, color=:orange, label="", foreground_color_legend=nothing);

expected_bites_treated = plot(a₀_bites, bites_age_treated_ageind, ribbon=bites_age_treated_ageind_unc, fillalpha=0.1, line=(3, :dot), color=:purple4, legendtitle="Treated", legend=:topright, label="", foreground_color_legend=nothing, xlabel=L"a_0~\textrm{(age infectious blood-meal taken in days)}", ylabel=L"E(z~|~\textrm{inf. blood-meal at}~a_0)");
plot!(a₀_bites, bites_age_treated_logistic, ribbon=bites_age_treated_logistic_unc, fillalpha=0.1, line=(3, :dash), color=:teal, label="", foreground_color_legend=nothing);
plot!(a₀_bites, bites_age_treated_gompertz, ribbon=bites_age_treated_gompertz_unc, fillalpha=0.1, line=3, color=:orange, label="", foreground_color_legend=nothing);

expected_bites_plots = [expected_bites_control, expected_bites_treated]

titleplotbites = plot([], grid=false, showaxis=false, leg=false, title=L"\textrm{Expected number of bites given infectious blood-meal at}~a_0", titlefont=20);

legendplot = plot([], label="Gompertz", color=:orange, grid=false, showaxis=false, line=3, legend=:outertop, legendfontsize=10);
plot!([], label="Age-Independent", color=:purple4, line=(3, :dot));
plot!([], label="Logistic", color=:teal, line=(3, :dash), foreground_color_legend=nothing);

l = @layout[a{0.001h}
            grid(1, 2)
            b{0.1h}]
step3 = plot(titleplotbites, expected_bites_plots..., legendplot, layout=l, ylims=(0, 6.1), xticks=(0:5:35), xlims=(0, 35), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10, show=true);
display(step3)


bites_control = scatter([1], [expected_bites_ageind_control], xlims=(0,4), ylims=(0,6), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", show=true, ms=5, m=:rect, label="", c=:purple4, title="Control", size=(400,380), ylabel=("Expected number of bites"));
scatter!([2], [expected_bites_logistic_control], ms=5, m=:rect, label="", c=:teal);
scatter!([3], [expected_bites_gompertz_control], ms=5, m=:rect, label="", c=:orange, tickdirection=:out);
# display(bites_control)

bites_treated = scatter([1], [expected_bites_ageind_treated], xlims=(0,4), ylims=(0,6), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", show=true, ms=5, m=:rect, label="", c=:purple4, title="Treated", size=(400,380));
scatter!([2], [expected_bites_logistic_treated], ms=5, m=:rect, label="", c=:teal);
scatter!([3], [expected_bites_gompertz_treated], ms=5, m=:rect, label="", c=:orange, tickdirection=:out);
# display(bites_treated)

bites_plots = [bites_control, bites_treated]

step4 = plot(bites_plots..., layout=grid(1,2), ylims=(0, 6), tickdirection=:out, fontfamily="Times", size=(1000, 500), tickfontsize=10, yguidefont=15, xguidefont=10, show=true);
display(step4)


# Violin plot
μ_control_ageind_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
    μ_control_ageind_errors[i] = μ_ageind_control[1] .± margin_error(fit_control_ageind, 0.01*i)[1]
end

μ_treated_ageind_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
    μ_treated_ageind_errors[i] = μ_ageind_treated[1] .± margin_error(fit_treated_ageind, 0.01*i)[1]
end


μ_control_logistic_errors1 = Array{Measurement}(undef, 100)
μ_control_logistic_errors2 = Array{Measurement}(undef, 100)
μ_control_logistic_errors3 = Array{Measurement}(undef, 100)
for i in collect(1:100)
      μ_control_logistic_errors1[i] = μ_logistic_control[1] .± margin_error(fit_control_logistic, 0.01*i)[1]
      μ_control_logistic_errors2[i] = μ_logistic_control[2] .± margin_error(fit_control_logistic, 0.01*i)[2]
      μ_control_logistic_errors3[i] = μ_logistic_control[3] .± margin_error(fit_control_logistic, 0.01*i)[3]
end

μ_treated_logistic_errors1 = Array{Measurement}(undef, 100)
μ_treated_logistic_errors2 = Array{Measurement}(undef, 100)
μ_treated_logistic_errors3 = Array{Measurement}(undef, 100)
for i in collect(1:100)
      μ_treated_logistic_errors1[i] = μ_logistic_treated[1] .± margin_error(fit_treated_logistic, 0.01*i)[1]
      μ_treated_logistic_errors2[i] = μ_logistic_treated[2] .± margin_error(fit_treated_logistic, 0.01*i)[2]
      μ_treated_logistic_errors3[i] = μ_logistic_treated[3] .± margin_error(fit_treated_logistic, 0.01*i)[3]
end


μ_control_gompertz_errors1 = Array{Measurement}(undef, 100)
μ_control_gompertz_errors2 = Array{Measurement}(undef, 100)
for i in collect(1:100)
      μ_control_gompertz_errors1[i] = μ_gompertz_control[1] .± margin_error(fit_control_gompertz, 0.01*i)[1]
      μ_control_gompertz_errors2[i] = μ_gompertz_control[2] .± margin_error(fit_control_gompertz, 0.01*i)[2]
end

μ_treated_gompertz_errors1 = Array{Measurement}(undef, 100)
μ_treated_gompertz_errors2 = Array{Measurement}(undef, 100)
for i in collect(1:100)
      μ_treated_gompertz_errors1[i] = μ_gompertz_treated[1] .± margin_error(fit_treated_gompertz, 0.01*i)[1]
      μ_treated_gompertz_errors2[i] = μ_gompertz_treated[2] .± margin_error(fit_treated_gompertz, 0.01*i)[2]
end


expected_bites_ageind_control_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
      expected_bites_ageind_control_errors[i] = expected_bites(μ_control_ageind_errors[i], "Age-Independent", "Erlang")
end

expected_bites_ageind_treated_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
      expected_bites_ageind_treated_errors[i] = expected_bites(μ_treated_ageind_errors[i], "Age-Independent", "Erlang")
end

expected_bites_logistic_control_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
      expected_bites_logistic_control_errors[i] = expected_bites([μ_control_logistic_errors1[i], μ_control_logistic_errors2[i], μ_control_logistic_errors3[i]], "Logistic", "Erlang")
end

expected_bites_logistic_treated_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
      expected_bites_logistic_treated_errors[i] = expected_bites([μ_treated_logistic_errors1[i], μ_treated_logistic_errors2[i], μ_treated_logistic_errors3[i]], "Logistic", "Erlang")
end

expected_bites_gompertz_control_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
      expected_bites_gompertz_control_errors[i] = expected_bites([μ_control_gompertz_errors1[i], μ_control_gompertz_errors2[i]], "Gompertz", "Erlang")
end

expected_bites_gompertz_treated_errors = Array{Measurement}(undef, 100)
for i in collect(1:100)
    expected_bites_gompertz_treated_errors[i] = expected_bites([μ_treated_gompertz_errors1[i], μ_treated_gompertz_errors2[i]], "Gompertz", "Erlang")
end

# Age-Independent
vc_control_ageind_errors = VectorialCapacity(expected_bites_ageind_control_errors)

vc_treated_ageind_errors = VectorialCapacity(expected_bites_ageind_treated_errors)

# Logistic
vc_control_logistic_errors = VectorialCapacity(expected_bites_logistic_control_errors)

vc_treated_logistic_errors = VectorialCapacity(expected_bites_logistic_treated_errors)

# Gompertz
vc_control_gompertz_errors = VectorialCapacity(expected_bites_gompertz_control_errors)

vc_treated_gompertz_errors = VectorialCapacity(expected_bites_gompertz_treated_errors)

# Differences
ageind_difference = VCCalculations.TreatmentDecrease(vc_control_ageind_errors, vc_treated_ageind_errors)
ageind_value = Measurements.value(ageind_difference[1])
ageind_unc = Measurements.uncertainty.(ageind_difference)
ageind_errors_a = ageind_value .+ ageind_unc
ageind_errors_b = reverse(ageind_value .- ageind_unc)[2:end]
ageind_errors = vcat(ageind_errors_a, ageind_errors_b)

logistic_difference = VCCalculations.TreatmentDecrease(vc_control_logistic_errors, vc_treated_logistic_errors)
logistic_value = Measurements.value(logistic_difference[1])
logistic_unc = Measurements.uncertainty.(logistic_difference)
logistic_errors_a = logistic_value .+ logistic_unc
logistic_errors_b = reverse(logistic_value .- logistic_unc)[2:end]
logistic_errors = vcat(logistic_errors_a, logistic_errors_b)

gompertz_difference = VCCalculations.TreatmentDecrease(vc_control_gompertz_errors, vc_treated_gompertz_errors)
gompertz_value = Measurements.value(gompertz_difference[1])
gompertz_unc = Measurements.uncertainty.(gompertz_difference)
gompertz_errors_a = gompertz_value .+ gompertz_unc
gompertz_errors_b = reverse(gompertz_value .- gompertz_unc)[2:end]
gompertz_errors = vcat(gompertz_errors_a, gompertz_errors_b)

# Plotting the violin plot
using StatsPlots

violinplot = violin(repeat([1], 199), ageind_errors, label="", alpha=0.3, tickdirection=:out, c=:purple4, ylims=(0, 100), xticks=((0:4), ["", "Age-Independent", "Logistic", "Gompertz", ""]), fontfamily="Times", ylabel="Percentage relative to control");
scatter!([1], [ageinderl], ms=5, m=:rect, label="", c=:purple4);
violin!(repeat([2], 199), logistic_errors, label="", alpha=0.3, c=:teal);
scatter!([2], [logisticerl], ms=5, m=:rect, label="", c=:teal);
violin!(repeat([3], 199), gompertz_errors, label="", alpha=0.2, c=:orange);
scatter!([3], [gompertzerl], ms=5, m=:rect, label="", c=:orange);
yticks!(0:20:100);
display(violinplot)
