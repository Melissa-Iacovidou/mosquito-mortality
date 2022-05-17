using XLSX, Plots, DataFrames

# Edit file paths accordingly

# Figures for initial data in the paper and cumulative mortality in the Supporting information section

df_control =  DataFrame(XLSX.readtable("../../Supporting_information/S1_Table.xlsx", "CumulMortalityControl")...);

df_treated =  DataFrame(XLSX.readtable("../../Supporting_information/S2_Table.xlsx", "CumulMortalityTreated")...);

df_survival_control = DataFrame(XLSX.readtable("../../Supporting_information/S5_Table.xlsx", "SurvivalControl")...);

df_survival_treated = DataFrame(XLSX.readtable("../../Supporting_information/S6_Table.xlsx", "SurvivalTreated")...);

# Extract data from dataframes
total_control = df_survival_control.Total;
total_treated = df_survival_treated.Total;

ages_control = df_control.Age;
ages_treated = df_treated.Age;

# Extract all the relevant data without the "Missing" values for each replicate
control_replicate_1 = 50 .- collect(skipmissing(df_control.Replicate1));
control_replicate_2 = 50 .- collect(skipmissing(df_control.Replicate2));
control_replicate_3 = 50 .- collect(skipmissing(df_control.Replicate3));
control_replicate_4 = 50 .- collect(skipmissing(df_control.Replicate4));
control_replicate_5 = 50 .- collect(skipmissing(df_control.Replicate5));

treated_replicate_1 = 50 .- collect(skipmissing(df_treated.Replicate1));
treated_replicate_2 = 50 .- collect(skipmissing(df_treated.Replicate2));
treated_replicate_3 = 50 .- collect(skipmissing(df_treated.Replicate3));
treated_replicate_4 = 50 .- collect(skipmissing(df_treated.Replicate4));
treated_replicate_5 = 50 .- collect(skipmissing(df_treated.Replicate5));

# Plot for control
control = plot(ages_control[1:length(control_replicate_2)], control_replicate_2, line=3, tickdirection=:out, ylim=(0,50), ylabel="alive mosquitoes", xaxis="age (in days)", xlim=(4,64), xticks=4:10:64, legend=false, c=:green, xguidefont=20);
plot!(ages_control[1:length(control_replicate_1)], control_replicate_1, line=3, c=:blue);
plot!(ages_control[1:length(control_replicate_3)], control_replicate_3, line=3, c=:red);
plot!(ages_control[1:length(control_replicate_4)], control_replicate_4, line=3, c=:orange);
plot!(ages_control[1:length(control_replicate_5)], control_replicate_5, line=3, c=:purple);
display(control)

# Plot for treated
treated = plot(ages_treated[1:length(treated_replicate_2)], treated_replicate_2, line=3, tickdirection=:out, ylim=(0,50), ylabel="", xaxis="age (in days)", xlim=(4,64), xticks=4:10:64, legend=false, c=:green, ytickfontcolor=:white, xguidefont=20);
plot!(ages_treated[1:length(treated_replicate_1)], treated_replicate_1, line=3, c=:blue);
plot!(ages_treated[1:length(treated_replicate_3)], treated_replicate_3, line=3, c=:red);
plot!(ages_treated[1:length(treated_replicate_4)], treated_replicate_4, line=3, c=:orange);
plot!(ages_treated[1:length(treated_replicate_5)], treated_replicate_5, line=3, c=:purple);
display(treated)

# Plot for proportion - control
prop_control = scatter(ages_control[1:length(control_replicate_2)], control_replicate_2 ./ 50, ms=5, tickdirection=:out, ylim=(-0.01,1.05), ylabel="proportion of alive mosquitoes", xaxis="age (in days)", xlim=(3.5,64.5), xticks=4:10:64, legend=false, c=:green, xguidefont=20, alpha=0.5);
scatter!(ages_control[1:length(control_replicate_3)], control_replicate_3 ./ 50, ms=5, c=:red, alpha=0.5);
scatter!(ages_control[1:length(control_replicate_4)], control_replicate_4 ./ 50, ms=5, c=:orange, alpha=0.5);
scatter!(ages_control[1:length(control_replicate_5)], control_replicate_5 ./ 50, ms=5, c=:purple, alpha=0.5);
plot!(ages_control, total_control ./ 200, line=5, legend=false, c=:dodgerblue);
display(prop_control)

# Plot for proportion - treated
prop_treated = scatter(ages_treated[1:length(treated_replicate_2)], treated_replicate_2 ./ 50, ms=5, tickdirection=:out, ylim=(-0.01,1.05), ylabel="", xaxis="age (in days)", xlim=(3.5,64.5), xticks=4:10:64, legend=false, c=:green, ytickfontcolor=:white, xguidefont=20, alpha=0.5);
scatter!(ages_treated[1:length(treated_replicate_3)], treated_replicate_3 ./ 50, ms=5, c=:red, alpha=0.5);
scatter!(ages_treated[1:length(treated_replicate_4)], treated_replicate_4 ./ 50, ms=5, c=:orange, alpha=0.5);
scatter!(ages_treated[1:length(treated_replicate_5)], treated_replicate_5 ./ 50, ms=5, c=:purple, alpha=0.5);
plot!(ages_treated, total_treated ./ 200, line=5, legend=false, c=:chocolate4);
display(prop_treated)

# Paper figures
legendplotcumul = plot([], showaxis = false, grid = false, color=:green, label = "Replicate 2", foreground_color_legend=nothing, legend=:outertop, line=3, legendfontsize=10, ticks=false);
plot!([], showaxis = false, grid = false, color=:blue, line=3, label="Replicate 1", ticks=false);
plot!([], color=:red, line=3, label="Replicate 3", ticks=false);
plot!([], color=:orange, line=3, label="Replicate 4", ticks=false);
plot!([], color=:purple, line=3, label="Replicate 5", ticks=false);

legendplotcumulprop = scatter([1], showaxis = false, grid = false, color=:green, label = "Replicate 2", foreground_color_legend=nothing, legend=:outertop, ms=5, legendfontsize=10, ticks=false);
scatter!([1], color=:red, ms=5, label="Replicate 3", ticks=false);
scatter!([1], color=:orange, ms=5, label="Replicate 4", ticks=false);
scatter!([1], color=:purple, ms=5, label="Replicate 5", ticks=false);
plot!([], color=:dodgerblue, line=3, label="Total Control", ticks=false);
plot!([], color=:chocolate4, line=3, label="Total Treated", ticks=false);


l = @layout [
    grid(1,2)
    a{0.2h}
];

cumulative_mortality = plot(control, treated, legendplotcumul, layout=l, title=["Control" "Treated" ""], tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true, size=(1800,900), titlefontsize=30);
display(cumulative_mortality)


l = @layout [
    grid(1,2)
    a{0.2h}
];

prop_cumulative_mortality = plot(prop_control, prop_treated, legendplotcumulprop, layout=l, title=["Control" "Treated" ""], tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", font="Times", show=true, size=(1800,900), titlefontsize=30);
display(prop_cumulative_mortality)


# Figure in Supporting information regarding feeding

df_feeding_control =  DataFrame(XLSX.readtable("../../Supporting_information/S3_Table.xlsx", "FeedingControl")...);

df_feeding_treated =  DataFrame(XLSX.readtable("../../Supporting_information/S4_Table.xlsx", "FeedingTreated")...);

ages_control = df_feeding_control.Age;
ages_treated = df_feeding_treated.Age;

# Extract all the relevant data without the "Missing" values for each replicate
feeding_control_replicate_1 = collect(skipmissing(df_feeding_control.Replicate1));
feeding_control_replicate_2 = collect(skipmissing(df_feeding_control.Replicate2));
feeding_control_replicate_3 = collect(skipmissing(df_feeding_control.Replicate3));
feeding_control_replicate_4 = collect(skipmissing(df_feeding_control.Replicate4));
feeding_control_replicate_5 = collect(skipmissing(df_feeding_control.Replicate5));

feeding_treated_replicate_1 = collect(skipmissing(df_feeding_treated.Replicate1));
feeding_treated_replicate_2 = collect(skipmissing(df_feeding_treated.Replicate2));
feeding_treated_replicate_3 = collect(skipmissing(df_feeding_treated.Replicate3));
feeding_treated_replicate_4 = collect(skipmissing(df_feeding_treated.Replicate4));
feeding_treated_replicate_5 = collect(skipmissing(df_feeding_treated.Replicate5));

# Plot for feeding - control
feeding_control = plot(ages_control[1:length(feeding_control_replicate_1)], feeding_control_replicate_1, line=3, tickdirection=:out, ylim=(0,1.0), ylabel="proportion fed from alive mosquitoes", xaxis="age (in days)", xlim=(4,60), xticks=4:8:60, legend=false, c=:blue, xguidefont=20);
plot!(ages_control[1:length(feeding_control_replicate_2)], feeding_control_replicate_2, line=3, c=:green);
plot!(ages_control[1:length(feeding_control_replicate_3)], feeding_control_replicate_3, line=3, c=:red);
plot!(ages_control[1:length(feeding_control_replicate_4)], feeding_control_replicate_4, line=3, c=:orange);
plot!(ages_control[1:length(feeding_control_replicate_5)], feeding_control_replicate_5, line=3, c=:purple);
display(feeding_control)

# Plot for feeding treated
feeding_treated = plot(ages_treated[1:length(feeding_treated_replicate_1)], feeding_treated_replicate_1, line=3, tickdirection=:out, ylim=(0,1), ylabel="", xaxis="age (in days)", xlim=(4,60), xticks=4:8:60, legend=false, c=:blue, ytickfontcolor=:white, xguidefont=20);
plot!(ages_treated[1:length(feeding_treated_replicate_2)], feeding_treated_replicate_2, line=3, c=:green);
plot!(ages_treated[1:length(feeding_treated_replicate_3)], feeding_treated_replicate_3, line=3, c=:red);
plot!(ages_treated[1:length(feeding_treated_replicate_4)], feeding_treated_replicate_4, line=3, c=:orange);
plot!(ages_treated[1:length(feeding_treated_replicate_5)], feeding_treated_replicate_5, line=3, c=:purple);
display(feeding_treated)

# Paper figure
legendplotfeed = plot([], showaxis = false, ticks=false, color=:blue, label = "Replicate 1", foreground_color_legend=nothing, legend=:outertop, line=3);
plot!([], color=:green, line=3, label="Replicate 2", legendfontsize=10);
plot!([], color=:red, line=3, label="Replicate 3");
plot!([], color=:orange, line=3, label="Replicate 4");
plot!([], color=:purple, line=3, label="Replicate 5");


l = @layout [
    grid(1,2)
    a{0.2h}
];

proportion_fed = plot(feeding_control, feeding_treated, legendplotfeed, layout=l, title=["Control" "Treated" ""], tick_direction=:out, yguidefont=25, xtickfont=20, ytickfont=20, tickfontfamily="Times", titlefontsize=30, font="Times", show=true, size=(1800,900));
display(proportion_fed)
