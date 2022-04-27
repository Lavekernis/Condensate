using StatsBase
using Base
using Statistics
using Plots
using ProgressMeter


#Simulation parameters
β = 0.2
N = 20
iterations = 10000
#-------------------------


function stateEnergy(state, energy_levels)
    return sum(state .* energy_levels)
end

function Step(state, states_numerator, energy_levels)
    """For given state returns next configuration and its energy"""
    accep = 1
    while true
        from_w_orbital = sample(states_numerator, Weights(state))
        inter_state = copy(state)
        inter_state[from_w_orbital] -= 1
        to_w_orbital = sample(states_numerator, Weights(inter_state .+ 1))
        final_state = copy(inter_state)
        final_state[to_w_orbital] += 1
        E_f, E_i = stateEnergy(final_state, energy_levels), stateEnergy(state, energy_levels)
        r = rand(Float32,1)
        if exp(-(E_f-E_i)*β) > r[1]
            return final_state, E_f, accep^(-1)
        end
        accep += 1
    end
end

function acceptance_plot()
    list_of_acceptance = []
    @showprogress 1 for cut_off in collect(1:1:30)
        states_vector = zeros(Int64, 1, 2*cut_off+1);
        states_vector[1,1] = N;
        accep_coef = []
        states_numerator = collect(1:1:2*cut_off+1)
        energy_levels = collect(-cut_off:1:cut_off).^2   
        for i in collect(1:1:iterations)
            temp = Step(states_vector[i,:], states_numerator, energy_levels)
            states_vector = [states_vector; temp[1]']
            accep_coef = [accep_coef; temp[3]]
        end
        list_of_acceptance = [list_of_acceptance; mean(accep_coef[500:end])]
    end
end