using StatsBase
using Base
using Statistics
using Plots
using ProgressMeter


#Simulation parameters
β = 0.2
cut_off = 10
N = 20
iterations = 1000
#-------------------------

#Vectors describing occupations
states_numerator = collect(1:1:2*cut_off+1)              # 1, 2, 3, ...
energy_levels = collect(-cut_off:1:cut_off).^2           # E(n_min), ..., 0, ..., E(n_max)
#----------------------------------------------

function stateEnergy(state)
    return sum(state .* energy_levels)
end

function Step(state, energy)
    """For given state returns next configuration and its energy"""
    accep = 1
    while true
        from_w_orbital = sample(states_numerator, Weights(state))
        inter_state = copy(state)
        inter_state[from_w_orbital] -= 1
        to_w_orbital = sample(states_numerator, Weights(inter_state .+ 1))
        final_state = copy(inter_state)
        final_state[to_w_orbital] += 1
        E_f, E_i = stateEnergy(final_state), energy
        r = rand(Float32,1)
        if exp(-(E_f-E_i)*β) > r[1]
            return final_state, E_f, accep^(-1)
        end
        accep += 1
    end
end


states_vector = zeros(Int64, 1, 2*cut_off+1);
states_vector[1,1] = N;
states_energy_vector = [stateEnergy(states_vector[1,:])]
accep_coef = []
@showprogress 1 for i in collect(1:1:iterations)
    local temp = Step(states_vector[i,:], states_energy_vector[i])
    global states_energy_vector = [states_energy_vector; temp[2]]
    global states_vector = [states_vector; temp[1]']
    global accep_coef = [accep_coef; temp[3]]
end

display(scatter(collect(1:1:iterations+1), states_energy_vector, 
plot_title = "N = $N, Cut-off = $cut_off, β = $β, iterations = $iterations", 
plot_titlefontsize = 13,
xaxis = "Iteration step",
yaxis = "⟨E⟩"))


energies_f_histogram = []
for (index,value) in enumerate(states_vector[iterations+1,:])
    for j in collect(1:1:value)
        global energies_f_histogram = [energies_f_histogram; energy_levels[index]]
    end
end

display(heatmap(collect(-cut_off:1:cut_off),collect(100:1:iterations),states_vector[100:iterations,:],
plot_title = "N = $N, Cut-off = $cut_off, β = $β, iterations = $iterations", 
plot_titlefontsize = 13,
xaxis = "Wave number",
yaxis = "Iteration"))

histogram(N .- states_vector[1000:end,cut_off + 1], 
bins = 0:1:N+1,
plot_title = "N = $N, Cut-off = $cut_off, β = $β, iterations = $iterations", 
plot_titlefontsize = 13,
xaxis = "x",
yaxis = "p(Nₑₓ = x) ", normalize = :probability, size = (1000, 1000))


# Analytic part ---------------------------------------------
inf = 1000

function En(j)
    return j^2
end

function P(N, Nₑₓ)
    A = 0
    B = 1
    C = 0
    for j in (1:1:inf)
        F = 0
        G = 1
        for k in (1:1:inf)
            if k ≠ j
                F += (1 - exp(-β*( En(j) - En(k) )))^(-1)
                G *= (1 - exp(-β*( En(k) - En(j) )))^(-2)
            end
        end
        A += exp(-β*En(j)*Nₑₓ)*(Nₑₓ + 1 + 2*F)*G
        B *= (1 - exp(-β*En(j)))^(-2)
        C += (exp(-β*En(j)*N)/(exp(β*En(j)) - 1))*(N + 1 + (1 - exp(-β*En(j)))^(-1) + 2*F)*G
    end
    return A/(B-C)
end

prob_list = []
@showprogress 1 for Nₑₓ in (0:1:N)
    global prob_list = [prob_list ;P(N, Nₑₓ)]
end
#------------------------------------------------

display(scatter!(collect(0:1:N), prob_list,  
plot_titlefontsize = 13,
xaxis = "Nₑₓ",
yaxis = "P(N, Nₑₓ)",
size = (1000, 1000)))

mean_ground_state_occ = mean(states_vector[500:end, cut_off + 1])
std_ground_state_occ = std(states_vector[500:end, cut_off + 1])
mean_accep_coef = mean(accep_coef[500:end])
println("⟨N₀⟩ = $mean_ground_state_occ, σₙ₀ = $std_ground_state_occ, p = $mean_accep_coef")
