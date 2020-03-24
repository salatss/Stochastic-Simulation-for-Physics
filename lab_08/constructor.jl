include("lattice_P.jl")
using Random

struct Ising{T<:AbstractFloat,D} 
    N::Int               # number of spin
    β::T                 # inverse temperature
    Λ::Lattice{D}        # Lattice
    spin::Vector{Int}    # spin config
    h::Vector{T}         # external field
    H::Vector{T}         # total field 
end

mutable struct Measures{T<:AbstractFloat}
    it::Int         # it contains the iteration step 
    β::T            # the inverse temperature
    ene::Vector{T}  # a vector of length nmeas containing the energy measured at time it 
    mag::Vector{T}  # a vector of length nmeas containing the magnetization measured at time it 
    conf::Vector{BitArray{1}} # a vector of length nmeas containing the spin configuration measured at time it 
end


spin2bool = x -> x > 0 #converte i -1 in false e 1 in true



#calcola l'energia
function energy(I::Lattice, h::Vector, spin::Vector)
    en=0
    close=I.neig
    for i in 1:(prod(I.dims))
            en=en-spin[i]*h[i]
        
        for j in close[i] 
            en=en-0.5*spin[i]*spin[j]
        end
    end
    return en
end

#è un extra, ha bisogno di avere meno argoemnti visto che genera lei stessa spin0 e h
function Ising_rand(I::Tuple, β)
    N=prod(I)
    Λ=Lattice(I)
    spin0=[rand([1,-1]) for i in 1:N]        # spin config
    h=   [rand()*10 for i in 1:N]          # external field
            
    H=zeros(N)
    temp=0
    for i in 1:N
        temp=0
        for j in Λ.neig[i]
           temp=temp+spin0[j]
        end
        H[i]=temp+h[i]
    end
    return Ising(N,β,Λ,spin0,h, H)
end

#come faccio ad avere T e D come nella struct?
function Ising(I::Tuple,h::Vector{Float64},β::Float64, spin0::Vector{Int})
    N=prod(I)
    Λ=Lattice(I)
    H=zeros(N)
    temp=0
    for i in 1:N
        temp=0
        for j in Λ.neig[i]
           temp=temp +spin0[j]
        end
        H[i]=temp+h[i]
    end
    return Ising(N,β,Λ,spin0,h, H)
end

#NON SERVE CALCOLARE TUTTA LA SOMMA CHE DA L'ENERGIA MA SOLO IL TERMINE CHE PROVIENE DALLO SPIN FLIPPATO
function onemcstep!(x::Ising, site::Int)
                             
    ΔE=2*x.H[site]*x.spin[site]

    r=rand()

    if exp(-x.β*ΔE)>r
        x.spin[site]*=-1

        close=x.Λ.neig[site]

        for i in close
            x.H[i]+=2*x.spin[site]
        end

    end
    
end

#aggiorna tutti con un ordine
function onemcsweep!(x::Ising)

    perm=randperm(x.N)
    
    for i in perm
        onemcstep!(x,i)
    end

end

function measures!(output::Measures , ising::Ising )
    output.β=ising.β
    output.ene[output.it]=sum( -(ising.H+ising.h).* 0.5 .* ising.spin )
    output.mag[output.it]=abs(sum(ising.spin)/ising.N)
    output.conf[output.it]=BitArray(map(spin2bool,ising.spin))

end



function mcising(I::NTuple{N,Int},β::T;  # linear size of the lattice + temperature
    h::Vector{T}=zeros(prod(I)), # external field
    x0::Vector{Int}=rand([-1,1],prod(I)), # initial configuration 
    nterm::Int=100,  # number of mcsweeps for initial thermalization
    nmeas::Int=1000, # number of measuen
    nsweep::Int=100, # number of mcsweeps between measurements
    ising=Ising(I,h,β,x0) ) where {T,N}

    output=Measures(0,β, zeros(nmeas), zeros(nmeas), Vector{BitArray{1}}(undef, nmeas ))

    for t in 1:nterm   # thermalization
        onemcsweep!(ising)
    end
    
    for n in 1:nmeas # loop for measurements   
        
        output.it+=1
        for t in 1:nsweep # number of mcsweeps between measurements
            onemcsweep!(ising)
        end
        measures!(output,ising)
    end

    return output
end
