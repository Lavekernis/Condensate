using StatsBase
using Base
using Statistics
using Plots
using ProgressMeter

termalization = 5000

function σNₑₓ(β, g)
    #Simulation parameters
    cut_off = 5
    N = 100
    iterations = 300000
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
            E_f, E_i = stateEnergy(final_state), energy
            r = rand(Float32,1)
            if exp(-(E_f-E_i)*β) > r[1]
                return final_state, E_f
            end
            accep += 1
        end
    end

    states_vector = [zeros(2*cut_off+1,2*cut_off+1)]
    states_vector[1][1,1] = N
    states_energy_vector = [stateEnergy(states_vector[1])]
    for i in collect(1:1:iterations)
        local temp = Step(states_vector[i], states_energy_vector[i])
        push!(states_energy_vector, temp[2])
        push!(states_vector, temp[1])
    end

    N₀_list = []
    for states in states_vector[termalization:end]
        push!(N₀_list, N - states[cut_off + 1, cut_off + 1])
    end
    return std(N₀_list)
end

function disp(g)
    σ_list = Vector{Float64}()
    T_list = collect(0.1:0.1:10)
    @showprogress 1 for T in T_list
        push!(σ_list, σNₑₓ(1/T, g))
    end
    return T_list, σ_list
end

for g in (0:0.01:0.01)
    T_list, σ_list = disp(g)
    display(scatter!(T_list, σ_list,
    xaxis = "T",
    yaxis = "σNₑₓ",
    size = (1000, 1000), label = "g = $g"))
end



# io = open("dane.txt", "w")
# for i in (1:1:length(σ_list))
#     write(io, string(T_list[i]), "\t", string(σ_list[i]), "\n")
# end
# close(io)

# scatter(T_list, σ_list,
# xaxis = "T",
# yaxis = "σNₑₓ",
# size = (1000, 1000))