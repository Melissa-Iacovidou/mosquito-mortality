using XLSX, Plots, DataFrames

# Edit file paths accordingly

# Figures for initial data in the paper and cumulative mortality in the Supporting information section

df_control =  DataFrame(XLSX.readtable("../Supporting_information/S1_Table.xlsx", "CumulMortalityControl")...);

df_treated =  DataFrame(XLSX.readtable("../Supporting_information/S2_Table.xlsx", "CumulMortalityTreated")...);

df_survival_control = DataFrame(XLSX.readtable("../Supporting_information/S5_Table.xlsx", "SurvivalControl")...);

df_survival_treated = DataFrame(XLSX.readtable("../Supporting_information/S6_Table.xlsx", "SurvivalTreated")...);

# Extract data from dataframes
total_control = 200 .- df_survival_control.Total;
total_treated = 200 .- df_survival_treated.Total;

ages_control = df_control.Age;
ages_treated = df_treated.Age;

# Extract all the relevant data without the "Missing" values for each replicate
cumul_mort_control_replicate_1 = df_control.Replicate1[1:18];
cumul_mort_control_replicate_2 = df_control.Replicate2[1:40];
cumul_mort_control_replicate_3 = df_control.Replicate3[1:40];
cumul_mort_control_replicate_4 = df_control.Replicate4[1:61];
cumul_mort_control_replicate_5 = df_control.Replicate5[1:37];

cumul_mort_treated_replicate_1 = df_treated.Replicate1[1:16];
cumul_mort_treated_replicate_2 = df_treated.Replicate2[1:25];
cumul_mort_treated_replicate_3 = df_treated.Replicate3[1:29];
cumul_mort_treated_replicate_4 = df_treated.Replicate4[1:24];
cumul_mort_treated_replicate_5 = df_treated.Replicate5[1:23];

# Plot for control
cumul_mort_control = plot(ages_control[1:40], cumul_mort_control_replicate_2, line=3, tickdirection=:out, ylim=(0,50), ylabel="cumulative mortality", xaxis="age (in days)", xlim=(4,64), xticks=4:10:64, legend=false, c=:green, xguidefont=20);
plot!(ages_control[1:18], cumul_mort_control_replicate_1, line=3, c=:blue);
plot!(ages_control[1:40], cumul_mort_control_replicate_3, line=3, c=:red);
plot!(ages_control[1:61], cumul_mort_control_replicate_4, line=3, c=:orange);
plot!(ages_control[1:37], cumul_mort_control_replicate_5, line=3, c=:purple);
# display(cumul_mort_control)

# Plot for treated
cumul_mort_treated = plot(ages_treated[1:25], cumul_mort_treated_replicate_2, line=3, tickdirection=:out, ylim=(0,50), ylabel="", xaxis="age (in days)", xlim=(4,64), xticks=4:10:64, legend=false, c=:green, ytickfontcolor=:white, xguidefont=20);
plot!(ages_treated[1:16], cumul_mort_treated_replicate_1, line=3, c=:blue);
plot!(ages_treated[1:29], cumul_mort_treated_replicate_3, line=3, c=:red);
plot!(ages_treated[1:24], cumul_mort_treated_replicate_4, line=3, c=:orange);
plot!(ages_treated[1:23], cumul_mort_treated_replicate_5, line=3, c=:purple);
# display(cumul_mort_treated)

# Plot for proportion - control
cumul_mort_prop_control = scatter(ages_control[1:40], cumul_mort_control_replicate_2 ./ 50, ms=5, tickdirection=:out, ylim=(-0.01,1.05), ylabel="proportion cumulative mortality", xaxis="age (in days)", xlim=(3.5,64.5), xticks=4:10:64, legend=false, c=:green, xguidefont=20, alpha=0.5);
scatter!(ages_control[1:40], cumul_mort_control_replicate_3 ./ 50, ms=5, c=:red, alpha=0.5);
scatter!(ages_control[1:61], cumul_mort_control_replicate_4 ./ 50, ms=5, c=:orange, alpha=0.5);
scatter!(ages_control[1:37], cumul_mort_control_replicate_5 ./ 50, ms=5, c=:purple, alpha=0.5);
plot!(ages_control, total_control ./ 200, line=5, legend=false, c=:dodgerblue);
# display(cumul_mort_prop_control)

# Plot for proportion - treated
cumul_mort_prop_treated = scatter(ages_treated[1:25], cumul_mort_treated_replicate_2 ./ 50, ms=5, tickdirection=:out, ylim=(-0.01,1.05), ylabel="", xaxis="age (in days)", xlim=(3.5,64.5), xticks=4:10:64, legend=false, c=:green, ytickfontcolor=:white, xguidefont=20, alpha=0.5);
scatter!(ages_treated[1:29], cumul_mort_treated_replicate_3 ./ 50, ms=5, c=:red, alpha=0.5);
scatter!(ages_treated[1:24], cumul_mort_treated_replicate_4 ./ 50, ms=5, c=:orange, alpha=0.5);
scatter!(ages_treated[1:23], cumul_mort_treated_replicate_5 ./ 50, ms=5, c=:purple, alpha=0.5);
plot!(ages_treated, total_treated ./ 200, line=5, legend=false, c=:chocolate4);
# display(cumul_mort_prop_treated)

# Figures for paper
controltitle = plot([], xlabel="Control", grid=false, axis=false, leg=false, xguidefont=30);
treatedtitle = plot([], xlabel="Treated", grid=false, axis=false, leg=false, xguidefont=30);

legendplotcumul = plot([], showaxis = false, grid = false, color=:green, label = "Replicate 2", foreground_color_legend=nothing, legend=:outertop, line=3, legendfontsize=10);
plot!([], showaxis = false, grid = false, color=:blue, line=3, label="Replicate 1");
plot!([], color=:red, line=3, label="Replicate 3");
plot!([], color=:orange, line=3, label="Replicate 4");
plot!([], color=:purple, line=3, label="Replicate 5");

legendplotcumulprop = scatter([], showaxis = false, grid = false, color=:green, label = "Replicate 2", foreground_color_legend=nothing, legend=:outertop, ms=5, legendfontsize=10);
scatter!([], color=:red, ms=5, label="Replicate 3");
scatter!([], color=:orange, ms=5, label="Replicate 4");
scatter!([], color=:purple, ms=5, label="Replicate 5");
plot!([], color=:dodgerblue, line=3, label="Total Control");
plot!([], color=:chocolate4, line=3, label="Total Treated");


l = @layout [
    grid(1,2){0.05h}
    grid(1,2)
    a{0.2h}
];

cumulative_mortality = plot(controltitle, treatedtitle, cumul_mort_control, cumul_mort_treated, legendplotcumul, layout=l, tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true, size=(1800,900));
display(cumulative_mortality)


l = @layout [
    grid(1,2){0.05h}
    grid(1,2)
    a{0.2h}
];

prop_cumulative_mortality = plot(controltitle, treatedtitle, cumul_mort_prop_control, cumul_mort_prop_treated, legendplotcumulprop, layout=l, tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true, size=(1800,900));
display(prop_cumulative_mortality)
# savefig("./Plots/PropCumulMortality.svg")


# Figure in Supporting information regarding feeding

df_feeding_control =  DataFrame(XLSX.readtable("../Supporting_information/S3_Table.xlsx", "FeedingControl")...);

df_feeding_treated =  DataFrame(XLSX.readtable("../Supporting_information/S4_Table.xlsx", "FeedingTreated")...);

ages_control = df_feeding_control.Age;
ages_treated = df_feeding_treated.Age;

# Extract all the relevant data without the "Missing" values for each replicate
feeding_control_replicate_1 = df_feeding_control.Replicate1[1:5];
feeding_control_replicate_2 = df_feeding_control.Replicate2[1:10];
feeding_control_replicate_3 = df_feeding_control.Replicate3[1:10];
feeding_control_replicate_4 = df_feeding_control.Replicate4[1:15];
feeding_control_replicate_5 = df_feeding_control.Replicate5[1:9];

feeding_treated_replicate_1 = df_feeding_treated.Replicate1[1:5];
feeding_treated_replicate_2 = df_feeding_treated.Replicate2[1:6];
feeding_treated_replicate_3 = df_feeding_treated.Replicate3[1:7];
feeding_treated_replicate_4 = df_feeding_treated.Replicate4[1:6];
feeding_treated_replicate_5 = df_feeding_treated.Replicate5[1:7];

# Plot for feeding - control
feeding_control = plot(ages_control[1:5], feeding_control_replicate_1, line=3, tickdirection=:out, ylim=(0,1.0), ylabel="proportion fed from alive mosquitoes", xaxis="age (in days)", xlim=(4,60), xticks=4:8:60, legend=false, c=:blue, xguidefont=20);
plot!(ages_control[1:10], feeding_control_replicate_2, line=3, c=:green);
plot!(ages_control[1:10], feeding_control_replicate_3, line=3, c=:red);
plot!(ages_control[1:15], feeding_control_replicate_4, line=3, c=:orange);
plot!(ages_control[1:9], feeding_control_replicate_5, line=3, c=:purple);
# display(feeding_control)

# Plot for feeding treated
feeding_treated = plot(ages_treated[1:5], feeding_treated_replicate_1, line=3, tickdirection=:out, ylim=(0,1), ylabel="", xaxis="age (in days)", xlim=(4,60), xticks=4:8:60, legend=false, c=:blue, ytickfontcolor=:white, xguidefont=20);
plot!(ages_treated[1:6], feeding_treated_replicate_2, line=3, c=:green);
plot!(ages_treated[1:7], feeding_treated_replicate_3, line=3, c=:red);
plot!(ages_treated[1:6], feeding_treated_replicate_4, line=3, c=:orange);
plot!(ages_treated[1:7], feeding_treated_replicate_5, line=3, c=:purple);
# display(feeding_treated)

# Figure for paper
controltitle = plot([], xlabel="Control", grid=false, axis=false, leg=false, xguidefont=30);
treatedtitle = plot([], xlabel="Treated", grid=false, axis=false, leg=false, xguidefont=30);

legendplotfeed = plot([], showaxis = false, grid = false, color=:blue, label = "Replicate 1", foreground_color_legend=nothing, legend=:outertop, line=3);
plot!([], showaxis = false, grid = false, color=:green, line=3, label="Replicate 2", legendfontsize=10);
plot!([], color=:red, line=3, label="Replicate 3");
plot!([], color=:orange, line=3, label="Replicate 4");
plot!([], color=:purple, line=3, label="Replicate 5");


l = @layout [
    grid(1,2){0.05h}
    grid(1,2)
    a{0.2h}
];

proportion_fed = plot(controltitle, treatedtitle, feeding_control, feeding_treated, legendplotfeed, layout=l, tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true, size=(1800,900));
display(proportion_fed)
