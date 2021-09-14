module VCCalculations

export surviving_EIP, no_bites, bites_age, expected_bites, VectorialCapacity, TreatmentDecrease

using QuadGK, Cuba, Measurements

σ = 0.09698706953718507 # incubation rate for EIP distribution - obtained from the ErlangEIP.nb mathematica file
α = 0.25 # bite rate
k = collect(0:50) # used for number of bites
λ = 31*σ # variable for the Erlang distribution

# Rethinking the Vectorial Capacity calculations

# Step 1: What is the probability of a mosquito surviving the EIP given it takes an infectious blood-meal at age a₀?

function surviving_EIP(μ, a₀, case, EIPdistribution, σ = σ, λ = λ)
    # setting up vector with results
    if typeof(μ) == Array{Float64,1}
        vector_surviving_EIP = Array{Float64,1}(undef, length(a₀))
    elseif typeof(μ) == Array{Measurement{Float64},1}
        vector_surviving_EIP = Array{Measurement{Float64},1}(undef, length(a₀))
    else
        return "Error: type of μ not supported"
    end

    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                vector_surviving_EIP = repeat(σ ./ (σ .+ μ), length(a₀))
            elseif EIPdistribution == "Erlang"
                vector_surviving_EIP = repeat((λ ./ (λ .+ μ)).^(λ/σ), length(a₀))
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                for (i, a₀) in enumerate(a₀)
                    vector_surviving_EIP[i] = quadgk(t -> σ * exp(-σ*t) * exp(-μ[1] * t) * ((1 + exp(μ[2] * (μ[3] - (a₀ + t)))) / (exp(μ[2] * (μ[3] - a₀)) + 1)) ^ (-μ[1]/μ[2]), 0, +Inf)[1]
                end
            elseif  EIPdistribution == "Erlang"
                for (i, a₀) in enumerate(a₀)
                    vector_surviving_EIP[i] = quadgk(t -> (λ^(λ/σ) * t^(λ/σ-1) * exp(-λ*t) / factorial((λ/σ)-1)) * exp(-μ[1] * t) * ((1 + exp(μ[2] * (μ[3] - (a₀ + t)))) / (exp(μ[2] * (μ[3] - a₀)) + 1)) ^ (-μ[1]/μ[2]), 0, +Inf)[1]
                end
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                for (i, a₀) in enumerate(a₀)
                    vector_surviving_EIP[i] = quadgk(t -> σ * exp(-σ*t) * exp(-((μ[1] * exp(a₀ * μ[2]) * (-1 + exp(μ[2] * t))) / μ[2])), 0, 10^3)[1]
                end
            elseif EIPdistribution == "Erlang"
                for (i, a₀) in enumerate(a₀)
                    vector_surviving_EIP[i] = quadgk(t -> (λ^(λ/σ) * t^(λ/σ-1) * exp(-λ*t) / factorial((λ/σ)-1)) * exp(-((μ[1] * exp(a₀ * μ[2]) * (-1 + exp(μ[2] * t))) / μ[2])), 0, 10^3)[1]
                end
            else
                return "Error: Distribution not supported"
            end
        end
    else
        return "Error: Invalid case"
    end
    return vector_surviving_EIP
end


# Step 2: How many bites will the mosquito take if it has survived the EIP?

struct NoBitesMeans{no_bites, means}
    bites::no_bites
    means::means
end

function no_bites(μ, a₁, case, α = α, k = k)
    # setting up matrix with results
    if typeof(μ) == Array{Float64,1}
        matrix_bites = Array{Float64}(undef, length(k), length(a₁))
        means = Array{Float64,1}(undef, length(a₁))
    elseif typeof(μ) == Array{Measurement{Float64},1}
        matrix_bites = Array{Measurement{Float64}}(undef, length(k), length(a₁))
        means = Array{Measurement{Float64},1}(undef, length(a₁))
    else
        return "Error: type of μ not supported"
    end

    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            for (i, a₁) in enumerate(a₁)
                matrix_bites[:, i] = μ .* α .^ k ./ (μ .+ α) .^ (k .+ 1)
            end
            means = repeat(α ./ μ, length(a₁))
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            for (i, k) in enumerate(k)
                for (j, a₁) in enumerate(a₁)
                    matrix_bites[i, j] = quadgk(a₂ -> (α ^ k * (a₂ - a₁) ^ k) / factorial(big(k)) * exp(-((α + μ[1]) * (a₂ - a₁))) * μ[1] / (1 + exp(μ[2] * (μ[3] - a₂))) * ((exp(μ[2] * (μ[3] - a₂)) + 1) / (exp(μ[2] * (μ[3] - a₁)) + 1)) ^ (-μ[1]/μ[2]), a₁, 10^5)[1]

                    means[j] = quadgk(a₂ -> -(((a₁ - a₂) * exp((μ[1] * ((a₁ - a₂) * μ[2] + log(1 + exp(μ[2] * (-a₁ + μ[3]))) - log(1 + exp(μ[2] * (-a₂ + μ[3])))))/μ[2]) * α * μ[1])/(1 + exp(μ[2] * (-a₂ + μ[3])))), a₁, 10^5)[1]
                end
            end
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            for (i, k) in enumerate(k)
                for (j, a₁) in enumerate(a₁)
                    matrix_bites[i, j] = quadgk(a₂ -> (α ^ k * (a₂ - a₁) ^ k) / factorial(big(k)) * μ[1] * exp(-(-a₁ + a₂)*α + a₂*μ[2] - (μ[1]*(-exp(a₁*μ[2]) + exp(a₂*μ[2])))/μ[2]), a₁, 10^3)[1]

                    means[j] = quadgk(a₂ -> -((a₁ - a₂) * α * μ[1] * exp((exp(a₁ * μ[2]) * μ[1] - exp(a₂ * μ[2]) *μ[1] + a₂ * μ[2]^2)/μ[2])), a₁, 10^3)[1]
                end
            end
        end
    else
        return "Error: Invalid case"
    end
    return NoBitesMeans(matrix_bites, means)
end


# Step 3: What is the expected number of bites a mosquito takes in its lifetime if it has taken an infectious blood-meal at age a₀?

function bites_age(μ, a₀, case, EIPdistribution, σ = σ, α = α, λ = λ)
    # setting up vector with results
    bites_age = Array{Any}(undef, length(a₀)) # "Any" because we might have measurements or only floats

    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                bites_age = repeat(σ ./ (σ .+ μ) .* α ./ μ, length(a₀))
            elseif EIPdistribution == "Erlang"
                bites_age = repeat(((λ ./ (λ .+ μ)).^(λ/σ)) .* α ./ μ, length(a₀))
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                for (i, a₀) in enumerate(a₀)
                    integrate(μ₁, μ₂, μ₃) = cuhre((a, f) -> f[1] = -((a[2] *exp(-((a[2] * μ₁ + a[1] * (μ₁ - 2 *a[2]*μ₁ + σ - a[2]*σ))/((-1 + a[1]) * (-1 + a[2])))) * (1 + exp(μ₂ * (-a₀ + μ₃)))^(μ₁/μ₂) * (1 + exp(μ₂ * (-a₀ + a[1]/(-1 + a[1]) + a[2]/(-1 + a[2]) + μ₃)))^(-((μ₁ + μ₂)/μ₂)) * α * μ₁ * σ /((-1 + a[1])^2 * (-1 + a[2])^3))), 2)[1][1]

                    if typeof(μ) == Array{Float64, 1}
                        result = integrate(μ[1], μ[2], μ[3])
                        bites_age[i] = result
                    elseif typeof(μ) == Array{Measurement{Float64},1}
                        result = @uncertain integrate(μ[1], μ[2], μ[3])
                        bites_age[i] = result
                    else
                        return "Error: Invalid type of μ"
                    end
                end
            elseif EIPdistribution == "Erlang"
                for (i, a₀) in enumerate(a₀)
                    integrate(μ₁, μ₂, μ₃) = cuhre((a, f) -> f[1] = ((a[1]/(1 - a[1]))^(λ/σ) * a[2] * exp((a[1] *(-1 + a[2]) * λ - a[2] * μ₁ + a[1] * (-1 + 2*a[2]) * μ₁)/((-1 + a[1]) * (-1 + a[2]))) * (1 + exp(μ₂ * (-a₀ + μ₃)))^(μ₁/μ₂) * (1 + exp(μ₂ * (-a₀ + a[1]/(-1 + a[1]) + a[2]/(-1 + a[2]) + μ₃)))^(-((μ₁ + μ₂)/μ₂)) * α * λ^(λ/σ) * μ₁)/((-1 + a[1]) * a[1] * (-1 + a[2])^3 * factorial(-1 + λ/σ)), 2)[1][1]

                    if typeof(μ) == Array{Float64, 1}
                        result = integrate(μ[1], μ[2], μ[3])
                        bites_age[i] = result
                    elseif typeof(μ) == Array{Measurement{Float64},1}
                        result = @uncertain integrate(μ[1], μ[2], μ[3])
                        bites_age[i] = result
                    else
                        return "Error: Invalid type of μ"
                    end
                end
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                for (i, a₀) in enumerate(a₀)
                    integrate(μ₁, μ₂) = cuhre((a, f) -> f[1] = -((a[2] * exp((exp(a₀*μ₂) * μ₁)/μ₂ - (exp((a₀ + a[1]/(1 - a[1]) + a[2]/(1 - a[2])) * μ₂) * μ₁)/μ₂ + (a₀ + a[1]/(1 - a[1]) + a[2]/(1 - a[2])) * μ₂ + (a[1] * σ)/(-1 + a[1])) * α * μ₁ * σ)/((-1 + a[1])^2 * (-1 + a[2])^3)), 2)[1][1]

                    if typeof(μ) == Array{Float64, 1}
                        result = integrate(μ[1], μ[2])
                        bites_age[i] = result
                    elseif typeof(μ) == Array{Measurement{Float64},1}
                        result = @uncertain integrate(μ[1], μ[2])
                        bites_age[i] = result
                    else
                        return "Error: Invalid type of μ"
                    end
                end
            elseif EIPdistribution == "Erlang"
                for (i, a₀) in enumerate(a₀)
                    integrate(μ₁, μ₂) = cuhre((a, f) -> f[1] = ((a[1]/(1 - a[1]))^(λ/σ) * a[2] * exp((a[1]*(λ - μ₂))/(-1 + a[1]) + (exp(a₀ * μ₂) * μ₁ - exp((a₀ + a[1]/(1 - a[1]) + a[2]/(1 - a[2])) * μ₂) * μ₁ + (a₀ + a[2]/(1 - a[2])) * μ₂^2)/ μ₂) * α * λ^(λ/σ) *μ₁)/((-1 + a[1]) * a[1] * (-1 + a[2])^3 * factorial(λ/σ - 1)), 2)[1][1]

                    if typeof(μ) == Array{Float64, 1}
                        result = integrate(μ[1], μ[2])
                        bites_age[i] = result
                    elseif typeof(μ) == Array{Measurement{Float64},1}
                        result = @uncertain integrate(μ[1], μ[2])
                        bites_age[i] = result
                    else
                        return "Error: Invalid type of μ"
                    end
                end
            else
                return "Error: Distribution not supported"
            end
        end
    else
        return "Error: Invalid case"
    end
    return bites_age
end


# Step 4: What is the expected number of (infectious) bites a mosquito takes in its lifetime?

function expected_bites(μ, case, EIPdistribution, σ = σ, α = α, λ = λ)
    if case == "Age-Independent"
        if length(μ) != 1
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                expected_bites = σ / (σ + μ[1]) * α / μ[1]
            elseif EIPdistribution == "Erlang"
                expected_bites = ((λ / (λ + μ[1]))^(λ/σ)) * α / μ[1]
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Logistic"
        if length(μ) != 3
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                integrate1(μ₁, μ₂, μ₃) = cuhre((a, f) -> f[1] = -((a[3] * exp((a[1] * α)/(-1 + a[1]) - (a[3]*μ₁ + a[2] * (μ₁ - 2*a[3]*μ₁ + σ - a[3]*σ))/((-1 + a[2]) * (-1 + a[3]))) * (1 + exp(μ₂ * (a[1]/(-1 + a[1]) + μ₃)))^(μ₁/μ₂) * (1 + exp(μ₂ * (a[1]/(-1 + a[1]) + a[2]/(-1 + a[2]) + a[3]/(-1 + a[3]) + μ₃)))^(-((μ₁ + μ₂)/μ₂)) * α^2 * μ₁ * σ)/((-1 + a[1])^2 * (-1 + a[2])^2 * (-1 + a[3])^3)), 3)[1][1]

                if typeof(μ) == Array{Float64, 1}
                    expected_bites = integrate1(μ[1], μ[2], μ[3])
                elseif typeof(μ) == Array{Measurement{Float64},1}
                    expected_bites = @uncertain integrate1(μ[1], μ[2], μ[3])
                else
                    return "Error: Invalid type of μ"
                end
            elseif EIPdistribution == "Erlang"
                integrate2(μ₁, μ₂, μ₃) = cuhre((a, f) -> f[1] = ((a[2]/(1 - a[2]))^(λ/σ) * a[3] * exp((a[3]*μ₁ + a[2]*(λ - a[3]*λ + μ₁ - 2*a[3]*μ₁) + a[1]*((-1 + a[2])*(-1 + a[3]) * α + a[2] * (-1 + a[3]) * λ -a[3]*μ₁ + a[2]*(-1 + 2*a[3]) * μ₁))/((-1 + a[1]) * (-1 + a[2]) * (-1 + a[3]))) * (1 + exp(μ₂ * (a[1]/(-1 + a[1]) + μ₃)))^(μ₁/μ₂) * (1 + exp(μ₂ * (a[1]/(-1 + a[1]) + a[2]/(-1 + a[2]) + a[3]/(-1 + a[3]) + μ₃)))^(-((μ₁ + μ₂)/μ₂)) * α^2* λ^(λ/σ) *μ₁)/((-1 + a[1])^2 * (-1 + a[2]) * a[2] * (-1 + a[3])^3 * factorial(-1 + λ/σ)), 3)[1][1]

                if typeof(μ) == Array{Float64, 1}
                    expected_bites = integrate2(μ[1], μ[2], μ[3])
                elseif typeof(μ) == Array{Measurement{Float64},1}
                    expected_bites = @uncertain integrate2(μ[1], μ[2], μ[3])
                else
                    return "Error: Invalid type of μ"
                end
            else
                return "Error: Distribution not supported"
            end
        end
    elseif case == "Gompertz"
        if length(μ) != 2
            return "Error: Invalid length of μ"
        else
            if EIPdistribution == "Exponential"
                integrate3(μ₁, μ₂) = cuhre((a, f) -> f[1] = (-μ₁ * a[3] * α^2 * exp((a[1] *(α - μ₂))/(-1 + a[1]) - (a[3] * μ₂)/(-1 + a[3]) - (μ₁ * exp(-((a[1]*μ₂)/(-1 + a[1]))) * (-1 + exp(((a[2] + a[3] - 2 * a[2] * a[3]) * μ₂) / ((-1 + a[2]) * (-1 + a[3])))))/μ₂ + (a[2] * (-μ₂ + σ))/(-1 + a[2])) * σ)/((-1 + a[1])^2 * (-1 + a[2])^2 * (-1 + a[3])^3), 3, 1, atol=1e-12, rtol=1e-10)[1][1]

                if typeof(μ) == Array{Float64, 1}
                    expected_bites = integrate3(μ[1], μ[2])
                elseif typeof(μ) == Array{Measurement{Float64},1}
                    expected_bites = @uncertain integrate3(μ[1], μ[2])
                else
                    return "Error: Invalid type of μ"
                end
            elseif EIPdistribution == "Erlang"
                integrate4(μ₁, μ₂) = cuhre((a, f) -> f[1] = ((a[2]/(1 - a[2]))^(λ/σ) * a[3] * exp((a[1] * (α - μ₂))/(-1 + a[1]) + (a[2] * (λ - μ₂))/(-1 + a[2]) - (exp(-((a[1]*μ₂)/(-1 + a[1]))) * (-1 + exp(((a[2] + a[3] - 2*a[2]*a[3]) * μ₂)/((-1 + a[2]) * (-1 + a[3])))) * μ₁)/μ₂ - (a[3] * μ₂)/(-1 + a[3])) * α^2 * λ^(λ/σ) * μ₁)/((-1 + a[1])^2 * (-1 + a[2]) * a[2] * (-1 + a[3])^3 * factorial(λ/σ - 1)), 3, 1, atol=1e-12, rtol=1e-10)[1][1]

                if typeof(μ) == Array{Float64, 1}
                    expected_bites = integrate4(μ[1], μ[2])
                elseif typeof(μ) == Array{Measurement{Float64},1}
                    expected_bites = @uncertain integrate4(μ[1], μ[2])
                else
                    return "Error: Invalid type of μ"
                end
            else
                return "Error: Distribution not supported"
            end
        end
    else
        return "Error: Invalid case"
    end
    return expected_bites
end


# Vectorial Capacity and percentage decreased
# m, α, c, and b are constants so their values are not very important here as we check the TreatmentDecrease later, hence they cancel out. Have the function here for completeness.
function VectorialCapacity(exp_no_bites, m=50, α=0.25, c=0.91, b=0.24)
    return m * α * c * b .* exp_no_bites
end

function TreatmentDecrease(VC_control, VC_treated)
    reduction = VC_treated ./ VC_control .* 100
    return reduction
end

end
