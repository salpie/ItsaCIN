using StatsKit, Random, Distributions, Plots, StatsPlots, VegaLite, JLD2, Printf

### need to work on R and RMAX, check pop size, put in loop so can try different combinations of alpha, figure out how to plot tumour beast - maybe just save the cn and index of tumour beast instead of the whole array, and then calculate diversity inter and intra
##### structures # something means that when you get the cell divions there are x2 with ploidy 2.5 and x1 with 1.5... need to fix.

mutable struct Chrom
    CN::Array{Int64, 1} #8p gain and 13q gain
end

# mutable struct Chrfitness
#     optimum::Array{Int64, 1}
#     alpha::Array{Float64, 1}
# end

mutable struct CancerCell
    chrom::Chrom
    # chrfitness::Chrfitness
    N::Int64
    b::Float64
    d::Float64
    μ::Float64
    ploidy::Float64
end

function copychr(chromold)
    chromnew = chromold
    chromnew.CN = copy(chromold.CN)
    return chromnew
end

# function copychrfit(chromfitold)
#     chromfitnew = chromfitold
#     # chromfitnew.optimum = copy(chromfitold.optimum)
#     # chromfitnew.alpha = copy(chromfitold.alpha)
#     return chromfitnew
# end

function copycell(cancercellold::CancerCell)
  newcancercell::CancerCell = CancerCell(
  copychr(Chrom(cancercellold.chrom.CN)),
  # copychrfit(Chrfitness(cancercellold.chrfitness.optimum, cancercellold.chrfitness.alpha)),
  copy(cancercellold.N),
  copy(cancercellold.b),
  copy(cancercellold.d),
  copy(cancercellold.μ),
  copy(cancercellold.ploidy))
end

function wrongCN(cancercellold, Max_gain = 6, Max_loss = 1)
    verdict = 0
    if (maximum(cancercellold.chrom.CN) > Max_gain || minimum(cancercellold.chrom.CN) < Max_loss)
        verdict = 1
    else verdict = 0
    end
    return verdict
end

function halfweightedEuc(x,y,w)
    value = (x - y).^2 .* w
    return value
end

function optimumfitness(cancerCell, increasebirth)
    if increasebirth == true
            chr_weight = Float64[]
            for chr in 1:length(cancerCell.chrom.CN)
                push!(chr_weight, halfweightedEuc(optimumCN[chr], cancerCell.chrom.CN[chr], alpha[chr]))
            end

            dist = sqrt(sum(chr_weight)) #get weight eucledian

            smax = 0.3

            s = smax / (dist + 1)
            b = s + cancerCell.d
            d = b - s

            #d = cancercell.dinitial * (1 + propd * s)

            # b = (propb * b)  / (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            # d = (propd * d)  * (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            return b, d
    else
            chr_weight = Float64[]
            for chr in 1:length(cancerCell.chrom.CN)
                push!(chr_weight, halfweightedEuc(cancerCell.chrfitness.optimum[chr], cancerCell.chrom.CN[chr], cancerCell.chrfitness.alpha[chr]))
            end

            dist = sqrt(sum(chr_weight)) #get weight eucledian

            smax = cancerCell.b - cancerCell.d

            s = smax / (dist + 1)

            d = maximum([cancerCell.b - s, cancerCell.d])
            #d = cancercell.binitial - s
            b = s + d
            #d = cancercell.dinitial * (1 + propd * s)

            # b = (propb * b)  / (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            # d = (propd * d)  * (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            return b, d
    end
end

function updatePloidy(CancerCell)
    CancerCell.ploidy = mean(CancerCell.chrom.CN)
    return CancerCell
end

function adaptiveAlpha_optimum(CencerCell) # change the karyotype optimum and the alpha according to what is mutated

end

#    createalpha(alpha::Float64, N) = fill(alpha, N)
#    createalpha(alpha::Array{Float64, 1}, N) = alpha

function fillHeatmap(Tumour, time, tumour_beast, MaxPop)
    for cells in 1:length(Tumour)
#        tumour_beast[(MaxPop+1)-cells,time] = mean(Tumour[cells].chrom.CN)
        tumour_beast[cells,time] = mean(Tumour[cells].chrom.CN)
    end
return tumour_beast
end

function countmemb(itr)
    d = Dict{eltype(itr), Int}()
    for val in itr
        if isa(val, Number) && isnan(val)
            continue
        end
        d[val] = get!(d, val, 0) + 1
    end
    return d
end

function tt(births_array)
    r = similar(births_array)
    @inbounds for j = 1:size(births_array,2)
        @simd for i = 1:size(births_array,1)
            r[i,j] = births_array[1,j] * births_array[2,j] # fixed a typo here!
        end
    end
    r
end

function calcB(frame, val)
    g = similar(frame)
    @inbounds for j = 1:size(frame,2)
            g[2,j] = round(frame[2,j] * val) # fixed a typo here!
    end
    g
end

function allocateBirthDeaths(Tumour)
    cells = 3
    births = []
    number_births = round((length(Tumour)/100) * Tumour[cells].d)
    global number_births
    for i in 1:length(Tumour)
        append!(births, round(Tumour[i].b, sigdigits=2))
    end
    table = countmemb(births)
    births_array = hcat([[key, val] for (key, val) in table]...)
    frame = tt(births_array)
    frame[1,:] = births_array[1,:] #births on top row, number on bottom
    val = number_births/sum(frame[2,:])
    number_births_deaths = calcB(frame, val)
    number_births_deaths[1,:] = frame[1,:]
    population_number_stays_same = vcat(number_births_deaths, transpose(births_array[2,:]))

    return population_number_stays_same
end

function calculateStagnentPopulationBirths(Tumour, population_number_stays_same)
    deaths = unique(sort(sample(1:length(Tumour), Int(number_births), replace = false))) # delete this number of cells so we can have this number of births
    deleteat!(Tumour, unique(deaths))
    for i in 1:size(population_number_stays_same,2)
        while population_number_stays_same[2,i] > 0
            for s in 1:length(Tumour)
                roundedb  = round(Tumour[s].b, sigdigits=2)
                if roundedb == population_number_stays_same[1,i]
                    cancercellold = Tumour[s] # copy tumour cell
                    New_cancercell = copycell(cancercellold) # copy tumour cell
                    Tumour = push!(Tumour, New_cancercell)
                    population_number_stays_same[2,i] = population_number_stays_same[2,i] - 1 #1 less birth
                end
            end
        end
    end
    return Tumour
end
#
# anim = @animate for i ∈ 1:10
#     population_number_stays_same = allocateBirthDeaths(Tumour)
#     Tum = calculateStagnentPopulationBirths(Tumour, population_number_stays_same)
#     CNstates = hcat(map(x -> x.chrom.CN, Tum)...)'
#     heatmap(1:size(CNstates,2),
#         1:size(CNstates,1), CNstates,
#         c=:pu_or,
#         xlabel="chromosomes", ylabel="y values",
#         title="tumour cells", dpi=500)
# end



#death rate distribution
#μ = Poisson(0.02)
## add a mutation and see how this changes the dynamics
N = 22 # number of chromosomes
#optimumCN = [2,2,2,2,2,2,2,3,2,2,2,2,3,2,2,2,2,2,2,2,2,2]
optimumCN = [1,2,2,1,2,2,3,3,2,2,2,2,3,1,1,2,2,1,2,2,1,1]
global diploid = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
#alpha = [3,1,0.2,0.5,1,1,3,3,0.2,0.2,0.2,0.2,2,1,1,1,2,2,2,2,0.5,1]
alpha = [1,0.7,0.9,0.2,0.4,0.1,3,3,0.2,0.2,0.2,0.2,0.3,0.2,0.5,0.3,0.2,0.4,1,0.3,1,1]

μ_rate = Gamma(0.03)
#CN = fill(2, N)
CN = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
μ = rand(μ_rate, 1)
Max_gain = 6
Max_loss = 1
MaxPop = 10^5   # maximum population size of 1cm3 of epithelial cells, but double check 10^8
global alive = 0.00

Number_cells = []

### birth and death of normal diploid tumour
#birth rate distribution
brate = Normal(0.3,0.1)
b = rand(brate, 1)

#death rate distribution
drate = Normal(0.35,0.1)
d = rand(drate, 1)

#Initialise tumour
x = 1 # tumour size at first time point - one cancercell

Tumour = CancerCell[]
global Tumour
for i in 1:x
     push!(Tumour, CancerCell(Chrom(CN), 2, 1, 0.4, μ[1], 2.00))
#     push!(Tumour, CancerCell(Chrom(CN), Chrfitness(optimumCN, alpha), 2, b[1], d[1], μ[1], 2.00))
end

Time = 1500

# overall population
P = zeros(Int64,Time,1)

# keep a track of tumour dynamics
#tumour_beast = []
#tumour_beast = zeros(MaxPop+10^6, Time)
global list = []

global diploid_cells = zeros(Int64, Time)
global not_diploid_cells = zeros(Int64, Time)

# run the simulation, still need to add in function for Allee effecf
#function runSimulation(Time, tumour_beast, MaxPop, CN, Tumour, b, d)

function runSimulation(Time, MaxPop, CN, Tumour, b, d, diploid_cells, not_diploid_cells, optimumb, optimumd, Rmax)
    for time in 1:Time
        global number_births = 0
        global number_deaths = 0
        global number_nothing = 0
        global diploid_cell = 0
        global not_diploid_cell = 0
        global list = []

    println(time)
    println(length(Tumour))

    if (length(Tumour) <= MaxPop)
        total = 0
        for cell in 1:length(Tumour) # go through each cell in the tumour
            diploid_cell += Tumour[cell].chrom.CN == diploid # count the number of diploid cells
            not_diploid_cell += Tumour[cell].chrom.CN != diploid  # count the number of non-diploid cells

            counter = 1
            if cell == 1
                counter = 0
            end
#            if Tumour[cell].chrom.CN != diploid

#            end
#            r = rand(Uniform(rand_d, Rmax))
#            r = rand(Weibull(2.5,Rmax))
            r = rand(Uniform((Tumour[cell].d/5), Rmax))
#                        r = rand(Uniform(0, Rmax))
            # if time_counter < 3
            #     Rmax_first = Tumour[cell].b + Tumour[cell].d
            #     r = rand(Uniform(0.4, Rmax_first))
            # end
#            rb = rand(Uniform(0, optimumd))
            Tumour[cell].μ = rand(μ_rate, 1)[1]#assign mutation prob for each cell
            if time < 4
                Tumour[cell].μ = 1
            end
            global index = cell
            local cell
            if time >= 10
                increasebirth = true
                Tumour[cell].b = optimumfitness(Tumour[cell], increasebirth)[1]
                Tumour[cell].d = optimumfitness(Tumour[cell], increasebirth)[2]
            end
            if r < Tumour[cell].b

                # number_births = number_births+1 #increase number of births
                # if (Tumour[cell].chrom.CN == CN && cell > 3) #assign randomly from b and d distribution if diploid
                #     Tumour[cell].b = rand(b, 1)[1]
                #     Tumour[cell].d = rand(d, 1)[1]
                # end
                which_chrom = rand(1:N) #pick random chromosome to mutate
                if Tumour[cell].chrom.CN != diploid
                    if Tumour[cell].μ > 0.5 # if there is a mutation this cell division
                        Tumour[cell].chrom.CN[which_chrom]= Tumour[cell].chrom.CN[which_chrom]+1 # pick random chromosome to be gained
                    end
                end
                if Tumour[cell].chrom.CN == diploid
                    if Tumour[cell].μ > 0.7 # if there is a mutation this cell division
                        # which_chrom = rand(1:N) #pick random chromosome to mutate
                        Tumour[cell].chrom.CN[which_chrom]= Tumour[cell].chrom.CN[which_chrom]+1 # pick random chromosome to be gained
                    end
                end

                cancercellold = Tumour[cell] # copy tumour cell
                New_cancercell = copycell(cancercellold) # copy tumour cell

                if Tumour[cell].chrom.CN != diploid
                    if Tumour[cell].μ > 0.5 # if there is a mutation this cell division
                        New_cancercell.chrom.CN[which_chrom]= New_cancercell.chrom.CN[which_chrom]-2 # pick random chromosome to be gained
                    end
                end
                if Tumour[cell].chrom.CN == diploid
                    if Tumour[cell].μ > 0.7 # if there is a mutation this cell division
                        New_cancercell.chrom.CN[which_chrom]= New_cancercell.chrom.CN[which_chrom]-2 # pick random chromosome to be gained
                    end
                end

                deatholdcell = wrongCN(Tumour[cell]) # if the verdict is death by CN
                if (deatholdcell == 1)
                    # number_births = number_births-1 #increase number of deaths
                    # index = cell + number_births
                    # push!(list, index)
                end
                deathnewcell = wrongCN(New_cancercell)
                if (deathnewcell == 0)
                    counter = counter+1
                    number_births = number_births+1 #increase number of births
                    index_insert = cell + number_births
                    Tumour = insert!(Tumour, index_insert, New_cancercell) # stay alive and join the tumour
                end
            end
            total = total + counter
            if ( (Tumour[cell].b + Tumour[cell].d) <= r ) # if nothing -#                println("do nothing")                # do literally nothing
                number_nothing = number_nothing+1 #increase number of deaths
            end
            if (Tumour[cell].b <= r < Tumour[cell].b + Tumour[cell].d) # if b<r<b+d - then die #death event
#                number_deaths = number_deaths+1 #increase number of deaths
                index = total + 1
                push!(list, index)
            end
        end
    end
        push!(Number_cells, length(Tumour)) # count the number of tumour cells
        Tumour = deleteat!(Tumour, list)
        println(time)
        diploid_cells[time] = diploid_cell
        not_diploid_cells[time] = not_diploid_cell
end
    maxpop=true
    return Tumour, diploid_cells, not_diploid_cells, Number_cells
end

#sim = runSimulation(Time, tumour_beast, MaxPop, CN, Tumour, b, d)
for i in 1:500
    Tumour = CancerCell[]
    global Tumour
    CN = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
    for i in 1:x
        push!(Tumour, CancerCell(Chrom(CN), 2, 2, 0.4, 1, 2.00))
    end

    Tumour[1].chrom.CN = optimumCN
    increasebirth = true
    global optimumb = optimumfitness(Tumour[1], increasebirth)[1]
    global optimumd = optimumfitness(Tumour[1], increasebirth)[2]

    Tumour[1].chrom.CN = CN
                global Rmax = optimumb + optimumd
    sim = runSimulation(Time, MaxPop, CN, Tumour, b, d, diploid_cells, not_diploid_cells, optimumb, optimumd, Rmax)

    #@save "Tumour.jld2" Tumour
    #@save "diploid_cells.jld2" diploid_cells
    #@save "not_diploid_cells.jld2" not_diploid_cells

    p = plot(diploid_cells, label = "Diploid Cells", grid=false)
    plot!(p, not_diploid_cells, label = "Non-diploid Cells", grid=false)
    println(savefig(@sprintf("pop_sizes_try%03.d.jl", i)))

    CNstates = hcat(map(x -> x.chrom.CN, Tumour)...)'
    #@save "CNstates_new.jld2" CNstates
    if !isempty(CNstates)
        heatmap(1:size(CNstates,2),
            1:size(CNstates,1), CNstates,
    #        c=:pu_or,
            c=cgrad([:blue, :white, :pink, :red]),
            xlabel="chromosomes", ylabel="tumour cells",
            title="tumour", dpi=300)
    println(savefig(@sprintf("tumour_evo_try%03.d.jl", i)))
    end
end
