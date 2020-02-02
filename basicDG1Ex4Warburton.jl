using LinearAlgebra
using Plots

#=

1D scheme, with with lagragian first order basis functions on evenly spaced mesh.
Interpolation points are equal to element boundaries.

=#

# initial condition
f(x) = 0
g(t) = 1

function update_flux!(flux::Array{Float64}, u::Array{Float64}, t::Float64, a::Float64, alpha::Float64)
    #=

    for the far left element, the flux on the left side is the boundary condition,
    and the flux on the right is only dependent on the last element.
    
    =#
    #@show alpha
    flux[1] = g(t)
    flux[length(flux)] = u[length(u)]
    #print("t: ", t, "\n")
    #print("\tu: ", u, "\n")
    
    for i=2:2:length(flux)-2
        
        avg_t = (u[i] + u[i-1])
      
        jump_t = (1 - alpha) * (u[i] - u[i-1])
        flux_p = (a/2) * (avg_t + jump_t)
        #print("\t flux at (", i, "): ", flux_p, "\n")
        flux[i+1] = -flux_p
        flux[i] = flux_p
          
    end
end
                 
function main()
    let

        #=

        Parameters

        =#
        #wave-speed
        a = 1.0
        #growth-speed?
        b = 0
        #flux parameter
        alpha = 1.0
        # Left side of the interval
        l = -1
        # Right side of the interval
        r = 1
        # Number of elements
        nk = 3
        # order of approximations
        np = 2
        # number of basis functions
        nb = nk + 1
        # Length of each element
        h = (r - l)/nk
        # elements as a range
        x = l:h:r

        
        # Time parameters
        # Total time
        T = 1
        # time step
        t_step = .02
        # time interval
        t_int = t_step:t_step:T
        
        #=

        Assembing matrices and vectors

        =#

        # M matrix
        M = zeros(np,np)
        M[1,1] = 2/3 * (h/2)
        M[2,1] = 1/3 * (h/2)
        M[1, 2] = 1/3 * (h/2)
        M[1, 1] = 2/3 * (h/2)

        # S matrix
        S = zeros(np,np)
        S[1,1] = 1/2 * (h/2)
        S[2,1] = -1/2 * (h/2)
        S[1, 2] = -1/2 * (h/2)
        S[2, 2] = 1/2 * (h/2)

        #=
        Algebra to to get spatial discretization into form to do time discretization:
        
              Mᵏ∂uᵏ/∂t = (auᵏ)* ψ(xₗᵏ) - (auᵏ)* ψ(xᵣᵏ) + Sᵏ(auᵏ)+ Mᵏ(buᵏ)
              ⟹
              ∂uᵏ/∂t = (Mᵏ)⁻¹p + (aSᵏ + bMᵏ)uᵏ
        =#
 

        M⁻¹ = inv(M)
        E = M⁻¹ * (a .* S .+ b .* M)
        # inital fluxes for each element: nk + 1
        flux = Float64[]
        u = Float64[]
        x_graph = Float64[]
        for i in x           
            flux = vcat(flux, [-f(i) , f(i)])
            u = vcat(u, [f(i), f(i)])
            x_graph = vcat(x_graph, [i,i])
        end

        x_graph = x_graph[2:end-1]
        flux = flux[2:end-1]
        u = u[2:end-1]
        pushfirst!(u, g(0))

        
        plt = plot(legend=false, ylims = [-2, 2])
        #@show u

        # looping in time
        for t in t_int
            #loop over each element to update
            for k = 1:nk
                e_index = (2*k) - 1
                
                u[1] = g(t)
                plot!(x[k:k+1], u[e_index : e_index+1])
                dudt = M⁻¹ * flux[e_index : e_index+1] + E * u[e_index : e_index+1]
                @show flux[e_index : e_index+1]
                #@show t_step * dudt
                u[e_index:e_index+1] = u[e_index:e_index+1] + t_step * dudt
               
                
            end
            
            update_flux!(flux, u, t, a, alpha)
            display(plt)
            sleep(.1)
            #@show u
            #@show flux
            #@show t      
        end
    end
end


main()
