using StatsBase
using Base
using Statistics
using Plots
using ProgressMeter


#Simulation parameters
β = 0.2
cut_off = 5
N = 100
iterations = 10000
g = 0.1
#-------------------------

#Vectors describing occupation
states_numerator = Array{Vector}(undef, (2*cut_off+1)^2)
index = 1
for i in (1:1:2*cut_off+1)
    for j in (1:1:2*cut_off+1)
        states_numerator[(2*cut_off+1)*(i-1)+(j)] = [i,j]
    end
end
    
energy_levels = Array{Int128,2}(undef, 2*cut_off+1, 2*cut_off+1)
for i in (1:1:2*cut_off+1)
    for j in (1:1:2*cut_off+1)
        energy_levels[i,j] = (i-cut_off-1)^2 + (j-cut_off-1)^2 
    end
end

function stateEnergy(state)
    inter = 0
    for (i, i_n) in enumerate(state)
        for (j, j_n) in enumerate(state)
            if (i != j)
                inter += 1/2*g*i_n*j_n
            end
        end
    end 
    return sum(state .* energy_levels) + inter
end

function statemeanEnergy(state)
    return mean(state .* energy_levels)
end

function Step(state, energy)
    """For given state returns next configuration and its energy"""
    accep = 1
    while true
        test1 = vec(states_numerator)
        test2 = vec(state)
        from_w_orbital = sample(vec(states_numerator), Weights(vec(state)))
        inter_state = copy(state)
        inter_state[from_w_orbital[2], from_w_orbital[1]] -= 1
        to_w_orbital = sample(vec(states_numerator), Weights(vec(inter_state .+ 1)))
        final_state = copy(inter_state)
        final_state[to_w_orbital[2], to_w_orbital[1]] += 1
        E_f, E_i = stateEnergy(final_state), energy
        r = rand(Float32,1)
        if exp(-(E_f-E_i)*β) > r[1]
            return final_state, E_f, accep^(-1)
        end
        accep += 1
    end
end

states_vector = [zeros(2*cut_off+1,2*cut_off+1)]
states_vector[1][1,1] = N
states_energy_vector = [stateEnergy(states_vector[1])]
accep_coef = []
@showprogress 1 for i in collect(1:1:iterations)
    local temp = Step(states_vector[i], states_energy_vector[i])
    push!(states_energy_vector, temp[2])
    push!(states_vector, temp[1])
    push!(accep_coef, temp[3])
end

println(mean(accep_coef))
heatmap(states_vector[10001])

Nₑₓ_list = []
mean_energy = []
for state in states_vector[5000:end]
    push!(Nₑₓ_list, N - state[cut_off + 1, cut_off + 1])
    #push!(mean_energy, statemeanEnergy(state))
end

#display(scatter(mean_energy))

histogram(Nₑₓ_list, 
bins = 0:1:N+1,
plot_title = "N = $N, Cut-off = $cut_off, β = $β, iterations = $iterations", 
plot_titlefontsize = 13,
xaxis = "x",
yaxis = "p(Nₑₓ = x) ", normalize = :probability, size = (1000, 1000))

# function Z(N, β)
#     if N == 0
#         return 1
#     end
#     if N == 1
#         sum = 0
#         for e in energy_levels
#             sum += exp(-β*e)
#         end
#         return sum
#     end
#     sum = 0 
#     for n in (1:1:N)
#         sum += Z(1, n*β)*Z(N-n, β)
#     end
#     return 1/N*sum
# end

# expect = 0
# expect_square = 0
# for n in (1:1:N)
#     p_n = (Z(N-n, β)-Z(N-n-1, β))/Z(N,β)
#     global expect += p_n*n
#     global expect_square += p_n*n^2
# end
# print(expect_square - expect^2)

