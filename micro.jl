using StatsBase
using Base
using Statistics
using Plots
using ProgressMeter

termalization = 5000

#Simulation parameters
N = 100
cut_off = 5
iterations = 5000000
#-------------------------

function σNₑₓ(β, g)


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

    function stateEnergy(state, f_state, energy, from, to)
        ΔE = energy_levels[to[1], to[2]] - energy_levels[from[1], from[2]]
        Δinter = g/2*(state[from[2],from[1]]^2 + state[to[2],to[1]]^2 - f_state[from[2],from[1]]^2 - f_state[to[2],to[1]]^2)
        return energy + ΔE + Δinter
    end


    function Step(state, energy)
        """For given state returns next configuration and its energy"""
        accep = 1
        while true
            from_w_orbital = sample(vec(states_numerator), Weights(vec(state)))
            inter_state = copy(state)
            inter_state[from_w_orbital[2], from_w_orbital[1]] -= 1
            to_w_orbital = sample(vec(states_numerator), Weights(vec(inter_state .+ 1)))
            final_state = copy(inter_state)
            final_state[to_w_orbital[2], to_w_orbital[1]] += 1
            E_f, E_i = stateEnergy(state, final_state, energy, from_w_orbital, to_w_orbital), energy
            r = rand(Float32,1)
            if exp(-(E_f-E_i)*β) > r[1]
                return final_state, E_f
            end
            accep += 1
        end
    end

    states_vector = [zeros(2*cut_off+1,2*cut_off+1)]
    states_vector[1][cut_off + 1,cut_off + 1] = N
    states_energy_vector = [0.0]
    for i in collect(1:1:iterations)
        local temp = Step(states_vector[i], states_energy_vector[i])
        push!(states_energy_vector, temp[2])
        push!(states_vector, temp[1])
    end


    return states_vector, states_energy_vector
end

# function disp(g = 0.0)
#     T_list = collect(5:0.1:9)
#     states_matrix = Vector{Array}()
#     energy_matrix = Vector{Vector}()
#     @showprogress 1 for T_i in T_list
#         σ, state_vector, energy_vector = σNₑₓ(1/T_i, 0)
#         push!(states_matrix, state_vector)
#         push!(energy_matrix, energy_vector)
#     end
    
# end
state_vector, energy_vector = σNₑₓ(1/4, 0)
# display(histogram(energy_vector,
#     plot_titlefontsize = 13,
#     xaxis = ("E", 0:10:300),
#     yaxis = "N", size = (1000, 1000)))
   
    # Nₑₓ_list = []
    # for states in states_vector[termalization:end]
    #     push!(Nₑₓ_list, N - states[cut_off + 1, cut_off + 1])
    # end

hist = fit(Histogram, energy_vector)
E_max = argmax(hist.weights)*(maximum(energy_vector)/length(hist.weights))

σ_list = Vector{Float64}()
ΔE_list = (100:-1:0)
r_list = Vector{Int64}()
@showprogress 1 for ΔE in ΔE_list
    Nₑₓ_list = Vector{Int64}()
    r = 0
    for (i,E_i) in enumerate(energy_vector)
        if abs(E_i-E_max) < ΔE
            push!(Nₑₓ_list, N - state_vector[i][cut_off + 1, cut_off + 1])
            r += 1
        end
    end
    push!(σ_list, std(Nₑₓ_list))
    push!(r_list, r)
end
display(plot(ΔE_list, σ_list))
display(plot(ΔE_list, r_list))
