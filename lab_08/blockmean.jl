include("lattice_P.jl")
include("constructor.jl")
using Statistics


function blockmean(data::Vector{Measures{T}}; field::Symbol=:mag) where T<: AbstractFloat

    L=length(data);
    N=data[1].it; index=[];
    
    append!(index,N) 

    i=N
    while i>1   #creo il vettore delgi indici
        i/=2
        i=Int(ceil(i))
        append!(index,i)
    end
    index=reverse(index)
    
    MU=zeros(L,length(index)-1)

    for l in 1:L
        for i in 1:length(index)-1 #per non andare out of bounds visto che userò i+1
        #faccio la media dei dati compresi tra gli indici index[i] e index[i+1]
        #utilizzo getfield che mi permette di calcolare solo il campo specificato in input
            MU[l,i]=mean(getfield(data[l],field)[ index[i]:index[i+1]-1 ])  #il meno uno lo metto per non far ripetere nessun elemento
        end
    end

    t=index[2:end]
    μ=mean(MU,dims=1)
    σ=std(MU,dims=1)

    return t,μ',σ'
    
end

####################################################################
#AUTOCORRELAZIONE
####################################################################

function autocorrelation(data::Vector{Measures{T}},ΔT::Int; field::Symbol=:mag)  where T<:AbstractFloat

    M=length(data); n=length(getfield(data[1],field)); 
    μ=σ=zeros(ΔT+1) #perchè voglio partire da zero
    covs=zeros(ΔT+1,M)

    #nel caso di field=conf cè il problema che a ogni tempo ho un vettore e non uno scalare,
    #quindi devo creare un unico lungo vettore per calcolarne la covarianza
    if field==:conf
        L=length(data[1].conf[1]);
        uno=BitArray{1}(undef,L*n)  
        
        for Δt in 0:ΔT

            for i in 1:M
                
                for j in 1:n    #in questo ciclo faccio il collage di bitarray per ottenere un unico lungo vettore
                    uno[(j-1)*L+1 : j*L ]=data[i].conf[j] 
                end
                covs[Δt+1,i]=cor(uno[1:end-Δt*L], uno[1+Δt*L:end] ) #i am using pearson coefficient
            end
        end
    else
        for Δt in 0:ΔT

            for i in 1:M
                uno=getfield(data[i],field)
                covs[Δt+1,i]=cor(uno[1:end-Δt], uno[1+Δt:end] ) #i am using pearson coefficient
            end
        end
    end
    μ=mean(covs,dims=2)
    σ= std(covs,dims=2)
    return μ, σ
end