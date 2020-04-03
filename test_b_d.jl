using StatsKit, Random, Distributions, TimeSeries, Plots, StatsPlots, VegaLite

#birth rate distribution
b = Normal(0.5,0.2)

#death rate distribution
d = Normal(0.45,0.2)


d_death_probability = rand(d, 1000)
b_birth_probability = rand(b, 1000)

P = 10^6  # maximum population size

T = 2000
x = zeros(Int64,T,100)
x[1, :] .= 2.0
#x[:20]

for column in 1:100
    d_death_probability = rand(d, 3000)
    b_birth_probability = rand(b, 3000)
    for t in 1:(T-1)
        if 0 < x[t,column] < P-1
            # Is there a birth?
            birth = b_birth_probability[t]*x[t,column]
            # Is there a death?
            death = d_death_probability[t]*x[t,column]
            # We update the population size.
            x[t+1,column] = round(x[t,column] + ((1*birth) - (1*death)))
            # The evolution stops if we reach $0$ or $N$.
        else
            x[t+1,column] = x[t,column]
        end
    end
end

plot(x, title = "Birth Death Rates", xlabel="Time", ylabel="Population Size", grid=false, legend=false)
savefig("Birth_Death_detectable_size_2000_.pdf")

p = plot(b, label = "Birth Distribution", grid=false)
plot!(p, d, label = "Death Distribution", title="Birth and Death Distribution", grid=false)
savefig("Birth_Death_Distribution_detectable_size_2000_.pdf")

# calculate the number of deaths
number_deaths = x[10000,:]


#death rate distribution
#μ = Poisson(0.02)
μ = Gamma(0.02)
## add a mutation and see how this changes the dynamics

mutable struct CancerCell
    CN::Int64
    N::Int64
    b::Float64
    d::Float64
end

function copycell(cancercellold::cancercell)
    newcancercell::CancerCell = cancercell(copy(cancercellold.mutations), copy(cancercellold.death), copy(cancercellold.epnumber),copy(cancercellold.escaped))
end

#cancerCell = CancerCell(2, 1, 0.02, 0.01)
#Tumour=CancerCell[]
#Tumour[1] = CancerCell(2, 1, 0.02, 0.01)

x = 1 # I want the array to be of length 10
# create an uninitialized 1D array of size `x`
Tumour = Array{CancerCell}(undef, x) # or `Vector{Coords}(undef, x)`

# initialize the first tumour cancercell
for i = 1:x
    Tumour[i] = CancerCell(2, 1, 0.02, 0.01)
end

for time in 1:Time
    i = length(Tumour)
    for cell in 1:length(Tumour)
        if (Tumour[cell].b)
            cancercell_old = Tumour[cell]
            New_cancercell = copycell(cancercell_old)
            Tumour = push!(Tumour, New_cancercell)
        end
        if (Tumour[cell].d)
            Tumour = deleteat!(Tumour, cell)
        end
        if no change # do nothing
        end

    end
    P[time] = length(Tumour)
end

# copy functions - checkout?
push!(Tumour, CancerCell)
for cell in time
    for i in Tumour[cell,time]



########## new trial with a mutation added

time = 200

x = zeros(Int64,time,2)
#top row is no mutation
x[1, 1] = 2.0
#bottom row is mutation, mutation at time 300
x[70,2] = 2.0

P = 10^6  # maximum population size

### birth and death of normal tumour
#birth rate distribution
b = Normal(0.7,0.1)

#death rate distribution
d = Normal(0.5,0.1)

d_death_probability = rand(d, 1000)
b_birth_probability = rand(b, 1000)

### birth and death of tumour cells with copy number
#birth rate distribution
b_CN = Normal(0.79,0.1)

#death rate distribution
d_CN = Normal(0.5,0.1)

d_CNdeath_probability = rand(d_CN, 1000)
b_CNbirth_probability = rand(b_CN, 1000)

for i in 1:(time-1)
    if i >= 70
        if 0 < sum(x[i,1:2]) < P-1
            # Is there a birth?
            birth_CN = b_CNbirth_probability[i]*x[i,2]
            # Is there a death?
            death_CN = d_CNdeath_probability[i]*x[i,2]
            # We update the population size.
            x[i+1,2] = round(x[i,2] + ((1*birth_CN) - (1*death_CN)))
            # The evolution stops if we reach $0$ or $N$.
        else
            x[i+1,2] = x[i,2]
        end
    end
    if 0 < sum(x[i,1:2]) < P-1
        # Is there a birth?
        birth = b_birth_probability[i]*x[i,1]
        # Is there a death?
        death = d_death_probability[i]*x[i,1]
        # We update the population size.
        x[i+1,1] = round(x[i,1] + ((1*birth) - (1*death)))
        # The evolution stops if we reach $0$ or $N$.
    else
        x[i+1,1] = x[i,1]
    end
end

plot(x, title = "Birth Death Rates", xlabel="Time", ylabel="Population Size", grid=false)
savefig("Birth_Death_detectable_size_2000_withCN_high_birth_late_low_death.pdf")

p = plot(b, label = "Birth Distribution", grid=false)
plot!(p, d, label = "Death Distribution", title="Birth and Death Distribution", grid=false)
savefig("Birth_Death_Distribution_detectable_size_2000_withCN.pdf")


tumour_1 = [1:200]
tumour_2 = vcat(x[1:200,1],x[1:200,2])
diploid = ["chr13_diploid"]
gain = ["chr13_gain"]
g = vcat(repeat(diploid; outer=[200]), repeat(gain; outer=[200]))

a = areaplot(x = tumour_1, y = tumour_2, group = g, stacked = true)




#push a vector of Cells with each birth

#info for System - initialise the first CancerCell
cells = Cells(2, 1, 0, 0, 0)
T=3000
Total_Tumour = Array(CancerCells)

for column in 1:100
    d_death_probability = rand(d, 10000)
    b_birth_probability = rand(b, 10000)
    mutations_gain = rand(μ, 100000)
    mutations_loss = rand(μ, 100000)
    for t in 1:(T-1)
        if 0 < x[t,column] < P-1
            # Is there a birth?
            birth = b_birth_probability[t]*x[t,column]
            # Is there a death?
            death = d_death_probability[t]*x[t,column]
            # Is there a mutation?
            mutations_gain =
            mutations_loss = map(x -> rand(Poisson(x)), μ.loss) .> 0
            # We update the population size.
            x[t+1,column] = round(x[t,column] + ((1*birth) - (1*death)))
            # The evolution stops if we reach $0$ or $N$.
        else
            x[t+1,column] = x[t,column]
        end
    end
end


#function
CN=2
numberCellDivisions=1000
for i in 1:numberCellDivisions
    interger = rand(1:10)
    if interger == 4
        println("mutation")
    end
end


#function to add new mutations to cells based on μ
mutations_gain = map(x -> rand(Poisson(x)), μ.gain) .> 0
mutations_loss = map(x -> rand(Poisson(x)), μ.loss) .> 0

savefig("poisson_distribution.pdf")

mutable struct Chromosomes
    CN::Array{Int64, 1}
    N::Int64
    function Chromosomes(N, states = [])
        if isempty(states)
            CN = fill(2, N)
        else
            CN = states
        end
        return new(CN, N)
    end
endtypeof(2)
