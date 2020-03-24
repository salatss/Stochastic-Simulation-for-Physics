struct resND1
	x
	v
	t
end


function cromer_ND(f::Function, r0::Vector, v0::Vector, tlim::Vector, Δt::Number )

    m=1 #solve for unitary mass
	trange= tlim[1]:Δt:tlim[2]

	v=zeros( (length(trange), length(r0)) )
	r=zeros( (length(trange), length(r0)) )
   
	v[1,:]=v0
    r[1,:]=r0
   
	for i in 1:length(trange)-1
		v[i+1,:]=v[i,:]+f(r[i,:])/m*Δt
		r[i+1,:]=r[i,:]+v[i+1,:]*Δt

	end 
	
	res=resND1(r,v,trange)
	return res

end
