abstract type Solver end
struct Euler 		<: Solver end #Euler è tipo derivato di Solver
struct EulerCromer 	<: Solver end
struct MidPoint 	<: Solver end
struct LeapFrog 	<: Solver end

struct res
	x
	v
	t
end

function solve(solver::Type{Solver} ,f::Function,x0::Number,v0::Number,tlim::Vector,Δt::Number;verbose::Bool=false)
end

function solve(solver::Type{Euler} ,f::Function,x0::Number,v0::Number,tlim::Vector,Δt::Number;verbose::Bool=false)
	m=1 #solve for unitary mass
	trange= tlim[1]:Δt:tlim[2]
	v=Array{Float64}(undef,length(trange))
	x=Array{Float64}(undef,length(trange))
	v[1]=v0
	x[1]=x0

	for i in 1:length(trange)-1
		v[i+1]=v[i]+f(x[i])/m*Δt
		x[i+1]=x[i]+v[i]*Δt

	end 
	result=res(x, v,trange )
	
	return result
end

#function energy ()

function solve(solver::Type{EulerCromer},f::Function,x0::Number,v0::Number,tlim::Vector,Δt::Number;verbose::Bool=false)
	m=1 #solve for unitary mass
	trange= tlim[1]:Δt:tlim[2]
	v=Array{Float64}(undef,length(trange))
	x=Array{Float64}(undef,length(trange))
	v[1]=v0
	x[1]=x0

	for i in 1:length(trange)-1
		v[i+1]=v[i]+f(x[i])/m*Δt
		x[i+1]=x[i]+v[i+1]*Δt

	end 
	result=res(x, v,trange )
	
	return result
end

function solve(solver::Type{MidPoint},f::Function,x0::Number,v0::Number,tlim::Vector,Δt::Number;verbose::Bool=false)
	m=1 #solve for unitary mass
	trange= tlim[1]:Δt:tlim[2]
	v=Array{Float64}(undef,length(trange))
	x=Array{Float64}(undef,length(trange))
	v[1]=v0
	x[1]=x0

	for i in 1:length(trange)-1
		v[i+1]=v[i]+f(x[i])/m*Δt
		v_m=(v[i]+v[i+1])/2
		x[i+1]=x[i]+v_m*Δt
	end 
	result=res(x, v,trange )
	
	return result
end


function solve(solver::Type{LeapFrog},f::Function,x0::Number,v0::Number,tlim::Vector,Δt::Number;verbose::Bool=false)
	m=1 #solve for unitary mass
	trange= tlim[1]:Δt:tlim[2]
	x=zeros(length(trange))
	v=zeros(length(trange))

	vm=v0 +0.5*Δt*f(x0)
	v[1]=vm             #v[1/2]--->v[1] 
	x[1]=x0 + vm*Δt     #
	
	for i in 1:(length(trange)-1)
		v[i+1]=v[i] + Δt*f(x[i])
		x[i+1]=x[i] + Δt*(v[i+1])   #i è in verità semi intero
	end 
	
	
	result=res(x, v,trange )
	
	return result
end