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
