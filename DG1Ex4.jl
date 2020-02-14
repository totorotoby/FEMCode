using LinearAlgebra
using Plots

#=

1D scheme, with with lagragian first order basis functions on evenly spaced mesh.
Interpolation points are at the element boundaries.

=#

# initial condition

g(t) = 1


function update_flux!(flux::Array{Float64}, u::Array{Float64}, t::Float64, a::Float64, alpha::Float64)
    #=

    for the far left element, the flux on the left side is 0,
    and the flux on the right is only dependent on the last element.
    
    =#

    flux[1] = 1    
    flux[length(flux)] = u[length(u)]
    for i=2:length(flux)-1

        
        # getting the indexs for u⁺ and u⁻
        index_p = i*2 - 1
        index_m = index_p - 1

        avg_t = (u[index_p] + u[index_m])/2
        #@show u[index_p] , u[index_m]
        #@show avg_t
        
        jump_t = (1 - alpha)/2 * (u[index_p] - u[index_m])
        
        #@show jump_t
        flux[i] = a * (avg_t + jump_t)
        #@show i, flux[i]
       
    end
    #@show flux
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
        alpha = 0.0
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
        #@show h
        # elements as a range
        x = l:h:r
        
        #=
        Time parameters
        =#
        
        # Total time
        T = 5
        # time step
        t_step = .1
        # time interval
        t_int = t_step:t_step:T
        t_length = length(t_int)
        
        
        #=
        Assembing matrices and vectors
        =#
        
        # M matrix
        M = zeros(np,np)
        M[1,1] = 2/3 * (h/2)
        M[2,1] = 1/3 * (h/2)
        M[1, 2] = 1/3 * (h/2)
        M[2, 2] = 2/3 * (h/2)

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
              ∂uᵏ/∂t = (Mᵏ)⁻¹p + M⁻¹(aSᵏ + bMᵏ)uᵏ

        =#
        
        M⁻¹ = inv(M)
        E = M⁻¹ * (a .* S .+ b .* M)
        
        #@show M
        #@show M⁻¹
        #@show S
        #@show E
        #=
        inital u and x axis for graphing
        =#
        
        u = Float64[]
        x_graph = Float64[]
    
        for i in x
            u = vcat(u, [f(i), f(i)])
            x_graph = vcat(x_graph, [i,i])
        end
        
        x_graph = x_graph[2:end-1]
        u = u[2:end-1]
        

        #=
        inital flux
        =#

        flux = zeros(nb)
        update_flux!(flux, u, 0.0, a, alpha)


        @assert size(u,1) == 2*nk
        @assert size(x_graph,1) == 2*nk
        @assert size(flux,1) == nb
        #@assert flux[1] == 1
        @assert u[1] == f(0)
        
        #=
        main loops
        =#
        
        plt = plot(legend=false, ylims = [-2, 2])
        # looping in time
        for t in t_int
            print("\n")
            @show t
           
            #loop over each element to update
            for k = 1:nk

                # flux indices
                fl_index = k
                fr_index = fl_index + 1

                # solution indices
                ur_index = 2*k
                ul_index = ur_index - 1

                plot!(x_graph[ul_index : ur_index], u[ul_index : ur_index])
                
                # Debugging 
                #@show k
                #@show fl_index, fr_index
                #@show ul_index, ur_index
                
                dudt = M⁻¹ * [flux[fl_index], -flux[fr_index]] + E * u[ul_index : ur_index]
                if k == 1
                   
                    #@show [flux[fl_index],-flux[fr_index]]
                    #@show E
                    #@show u[ul_index : ur_index]
                    #@show ul_index, ur_index
                    #@show dudt
                    #@show E*u[ul_index : ur_index]
                    #@show M⁻¹*[flux[fl_index],-flux[fr_index]]
                    
                end
                u[ul_index : ur_index] = u[ul_index : ur_index] + t_step * dudt
     
            end
            
            update_flux!(flux, u, t, a, alpha)
            #@show u
            #@show flux
            #print("\n")
            display(plt)
            sleep(.2)
        end
    end
end


main()
