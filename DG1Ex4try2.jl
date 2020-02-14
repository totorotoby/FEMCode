using PyPlot

function update_flux!(flux::Array{Float64}, u::Array{Float64}, t::Float64, a::Float64, alpha::Float64)
    #=

    for the far left element, the flux on the left side is 0,
    and the flux on the right is only dependent on the last element.
    
    =#

    flux[1] = a*g(t)
    flux[length(flux)] = a * u[length(u)]
    for i=2:length(flux)-1

        
        # getting the indexs for u⁺ and u⁻
        index_p = i*2 - 1
        index_m = index_p - 1

        avg_t = (a * u[index_p] + a * u[index_m])/2
        #@show u[index_p] , u[index_m]
        #@show avg_t
        
        jump_t = a * (1 - alpha)/2 * (-u[index_p] + u[index_m])
       
        #@show jump_t
        flux[i] = (avg_t + jump_t)
        #@show i, flux[i]
       
    end
    #@show flux
end

function plotDG(data::Array{Array{Float64,1},1}, x_graph::Array{Float64, 1})

   
    fig = figure()
    ax = gca()
    ylims = (minimum([minimum(a) for a in data]), (maximum([maximum(a) for a in data])))
    axis(xlim=(x_graph[1], x_graph[end]), ylim = ylims)
    for t_step = 1:length(data)
        #@show data[t_step]
        line_segs = get_segs(x_graph, data[t_step])
        ls = matplotlib.collections.LineCollection(line_segs)
        ax.clear()
        ax.add_collection(ls)
        axis("image")
        sleep(.001)
        
    end  
end


function get_segs(x::Array{Float64,1}, y::Array{Float64, 1})
    
    segs = []
    pairs = collect(zip(x,y))
    
    for (i, pair) in enumerate(pairs)

        if mod(i,2) == 1
            push!(segs, [pairs[i], pairs[i+1]])
        end
        
    end
    
    return segs
end


# boundary condition
g(t) = 0#sin.(pi * t)
# inital condition
f(x) = 1#cos.(1/2*pi*x)

#=

nk = number of elements
t_step = length of timestep

=#

function assemble_global_matrix(m_inv, s, a, nk)

    A_loc = a.*m_inv*s

    fM_loc = [-a * m_inv[1,1] a*m_inv[1,1] ; -a * m_inv[2,1] a * m_inv[2,1]]
    @show fM_loc * [-1, 1]
    @show m_inv
    loc_size = size(m_inv, 1)
    
    A_g = zeros(2*nk, 2*nk)
        
    for i = 1:loc_size:size(A_g, 1)
            A_g[i:i+1,i:i+1] = A_loc
    end

    
    
    @show A_g
    
    
end


function advect(nk, t_step)
    
    let

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

        # solution vector:
        # length = 2*nk or np * nk            
        u = Float64[]
        for i in x
            u = vcat(u, [f(i),f(i)])
        end
        u = u[2 : end - 1]

        x_graph = Float64[]
        for i in x
            x_graph = vcat(x_graph, [i,i])
        end
        x_graph = x_graph[2 : end - 1]

        
        
        # inverse of mass matrix
        # dimensions : number of basis per element = np
        M⁻¹ = inv(M)
        #@show M⁻¹
        # "advection" matrix
        # dimensions : number of basis per element = np
        A = a .* M⁻¹ * S
        #@show A
        # flux vector:
        #  length = number of nodes = nb
        flux = zeros(nb)   
        
        for (t_i, t) in enumerate(t_int)
         
            #=
            for geting flux terms:
            k is the index of flux on the left of an element
            k+1 is the index of flux on the right of an element
            =#
            for k = 1:nk
                
                #indices for solution vector on element k
                ur_index = 2*k
                ul_index = ur_index - 1
                
                # the "flux term" for the left
                fl_term = -(a*u[ul_index] - flux[k])
                fr_term = (a*u[ur_index] - flux[k+1])

                
                dudt_s = - A * u[ul_index : ur_index] + b * u[ul_index : ur_index] + M⁻¹ * [fl_term, fr_term]
                #@show M⁻¹
                @show [fl_term, fr_term]
                @show  M⁻¹ * [fl_term, fr_term]
                @show u
                #@show  [fl_term, fr_term]
                #@show - A * u[ul_index : ur_index]
                print("\n")
                #dudt_w = A * u[ul_index : ur_index] + [flux[k] , -flux[k+1]]
                
                u[ul_index : ur_index] = u[ul_index : ur_index] + t_step * dudt_s
                
            end
            
            u_p = copy(u)
        
            push!(u_t, u_p)
            update_flux!(flux, u, t, a, alpha)
        
        end
        
        return u_t
       

    end
end



function convergence_table()

    
end



function main()
    let

           #=
        Parameters
        =#
        #wave-speed
        a = 1.0
        #growth-speed?
        b = 0.0
        #flux parameter
        alpha = 0.0
        # Left side of the interval
        l = -1
        # Right side of the interval
        r = 1
        # Number of elements
        #nk = 100
        # order of approximations
        np = 2
        # number of basis functions
        nb = nk + 1
        # Length of each element
        h = (r - l)/nk
        # elements as a range
        x = l:h:r
        
        #=
        Time parameters
        =#
        
        # Total time
        T = .02
        # time step
        #t_step = .1
        # time interval
        t_int = 0:t_step:T
        t_length = length(t_int)

        # solution to graph
        u_t = Array{Float64,1}[]

        
        
    end
end




nothing




