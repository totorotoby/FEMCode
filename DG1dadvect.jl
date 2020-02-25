using DelimitedFiles
using SparseArrays
using PyPlot
using Plots
using LinearAlgebra

function advec_global_strong_euler!(u::Array{Float64, 2}, G::SparseMatrixCSC{Float64, Int64},
                                    b::Array{Float64, 1}, t_int::StepRangeLen, step::Float64,
                                    g::Function, a::Float64)

    for (i, t) in enumerate(t_int)
        
        b[1] = a * g(t)
        dudt = G * u[i, :] + b
        u[i+1, :] = u[i, :] .+  dudt .* step
        
    end
end

function advec_local_strong_euler!(u::Array{Float64, 2}, A_l::Array{Float64, 2},
                                   F_l::Array{Float64, 2}, t_int::StepRangeLen,
                                   t_step::Float64, nk::Int64, g::Function, λ::Float64) 
    for (t_i, t) = enumerate(t_int)
        for k = 1:nk
            if k == 1
                dudt = A_l * u[t_i, 2*k - 1 : 2*k] + F_l * [g(t) , u[t_i, 1], u[t_i, 2], u[t_i, 3]]
            elseif k == nk
                dudt = A_l *  u[t_i, 2*k - 1 : 2*k] + F_l * [u[t_i, 2*k-2], u[t_i, 2*k-1], u[t_i, 2*k], u[t_i, 2*k]] 
            else                
                dudt = A_l * u[t_i, 2*k - 1 : 2*k] + F_l * u[t_i, 2*k-2 : 2*k+1]
            end
            u[t_i + 1, 2*k - 1 : 2*k] =  u[t_i, 2*k - 1 : 2*k] + t_step * dudt  
        end
    end
end

#=
function advec_global_rk!(u::Array{Float64, 2}, G::SparseMatrixCSC{Float64, Int64},
                          b::Array{Float64, 1}, t_int::StepRangeLen, step::Float64,
                          g::Function, a::Float64)
    
    for (i, t) in enmuerate(t_int)
        b[1] = a * g(t)
        k1 = G * u[i, :] + b
        k2 = 3
       
    end
end
=#

function get_local_mass!(M::Array{Float64, 2}, h::Float64)

    M[1,1] = 2/3 * (h/2)
    M[2,1] = 1/3 * (h/2)
    M[1, 2] = 1/3 * (h/2)
    M[2, 2] = 2/3 * (h/2)

end

function get_local_stiffness!(S::Array{Float64, 2}, h::Float64)
    
    S[1,1] = -1/2 #* (h/2)
    S[1,2] = 1/2 #* (h/2)
    S[2, 1] = -1/2 #* (h/2)
    S[2, 2] = 1/2 #* (h/2)
    
end

function get_local_flux!(F::Array{Float64, 2}, M_inv::Array{Float64, 2}, a::Float64, α::Float64)

    F_n = zeros(2,4)
    F_n[1,1] = a*α
    F_n[1,2] = -a*α
    F_n[2,3] = a*(1-α)
    F_n[2,4] = -a*(1-α)
    F[: , :] = M_inv*F_n
    
end

function get_global_advection!(A_g::SparseMatrixCSC{Float64, Int64}, A_l::Array{Float64, 2})
    
    for i=1 : 2 : size(A_g, 1)
        A_g[i :  i+1, i : i+1] = A_l
    end
end


function get_global_flux!(F_g::SparseMatrixCSC{Float64, Int64}, F_l::Array{Float64, 2})

    dim = size(F_g, 1)
    F_g[1 : 2 , 1 : 3] = F_l[ : , 2 : 4]
    F_g[dim-1 : dim, dim - 2 : dim] = F_l[ : , 1 : 3]

    for i = 3 : 2 : dim - 2
        F_g[i : i+1, i-1 : i+2] = F_l
    end

    
end
                                  

function get_inital_condition!(u::Array{Float64, 2}, f::Function , x::Array{Float64, 1}, nk::Int64)

    col_num = 2 * nk
    
    u[1,1] = f(x[1])
    
    u[1, col_num] = f(x[end])
    
    for (i, x_val) in enumerate(x[2 : end - 1])
        u[1, 2 * i] = f(x_val)
        u[1, (2 * i) + 1] = f(x_val)
    end
    
end


function write_plot(u::Array{Float64, 2}, x::Array{Float64, 1})
     
    open("num_data", "w") do io
        writedlm(io, u)
    end
    
    x_graph = Float64[]
    for i in x
        x_graph = vcat(x_graph, [i,i])
    end
    x_graph = x_graph[2 : end - 1]
    
    open("x_num_data", "w") do io
        writedlm(io, x_graph)
    end
    
end


function plot_DG(u_t::Array{Float64, 2}, x::Array{Float64, 1})

    x_graph = Float64[]
    for i in x
        x_graph = vcat(x_graph, [i,i])
    end
    x_graph = x_graph[2 : end - 1]

    fig = figure()
    ax = gca()
    ylims = (minimum(u_t), maximum(u_t))
    axis(xlim=(x_graph[1], x_graph[end]), ylim = ylims)
    for t_step = 1:4:size(u_t, 1)
        
        line_segs = get_segs(x_graph, u_t[t_step, :])
        ls = matplotlib.collections.LineCollection(line_segs)
        ax.clear()
        ax.add_collection(ls)
        axis("image")
        
        sleep(.0001)
        
    end
end


function write_analytic(x_an::Array{Float64,1}, t_int::Array{Float64,1}, a::Float64, λ::Float64)
   
    u_an =  Array{Float64, 2}(undef, size(t_int, 1), size(x_an, 1))
    u_a1!(u_an, a, λ, x_an, t_int)

    open("x_an_axis", "w") do io
        writedlm(io, x_an)
    end
    
    open("an_data", "w") do io
        writedlm(io, u_an)
    end
end


function get_segs(x_graph::Array{Float64,1}, u::Array{Float64, 1})
   
    segs = []
    pairs = collect(zip(x_graph,u))
    
    for (i, pair) in enumerate(pairs)

        if mod(i,2) == 1
            push!(segs, [pairs[i], pairs[i+1]])
        end       
    end
    return segs
end

function show_stability_region(G::Array{Float64,2}, h::Float64, t_step::Float64)

    E =eigen(t_step .* G)
    
    to_graph = E.values .* h
  
    plt = Plots.plot()
    scatter!(to_graph)
    return plt

end

function stability_test(G::Array{Float64,2}, h::Float64, t_step::Float64)

    E = eigen(t_step .* G).values
    sort!(E, by=abs)
    max_eig = E[end]
    @show abs(max_eig + 1)
    if (abs(max_eig + 1) < 1)
        return true
    else
        return false
    end
end


function u_a1!(u_a::Array{Float64,2}, a::Float64, λ::Float64, x::Array{Float64, 1},
               t::Array{Float64, 1})


    for (j, t_j) in enumerate(t)
        for (i, x_i) in enumerate(x)

            if (t_j*a < x_i + 1)
                u_a[j, i] = cos((π*(x_i -a*t_j))/2)*exp(λ*t_j)
            else
                u_a[j, i] = 0
            end
        end
    end
    
end



function convergence_test(a, Float64, λ::Float64, x::Array{Float64, 1},
                          t_int::Array{Float64, 1})

    #anayltic solution
    
    
end



function main()

    let
        ## Parameters ##
        
        # wave-speed
        a = 1.0
        # growth-speed?
        λ = 0.5
        # flux parameter
        α = 1.0
        # Left side of the interval
        l = -1
        # Right side of the interval
        r = 1
        # Number of elements
        nk = 100
        # order of approximations
        np = 2
        # Length of each element
        h = (r - l)/nk
        @show h
        @show a
        # elements as a range
        x = l:h:r
        # anayltic solution x values
        x_an = l:.01:r

        # filename for plotting
        filename = "plotting_data"
        
        ## Time parameters ##
        
        # Total time
        T = 1
        # time step
        t_step = .0001

        # time interval
        t_int = t_step:t_step:T
        t_num = length(t_int)
        
        # solution to graph
        u_t = Array{Float64,1}[]


        ## intial and boundary conditions ##

        f(x::Float64) = cos.((pi*x)/2)
        g(t::Float64) = 0

        ## Local Matrices ##

        # Mass matrix
        M = Array{Float64, 2}(undef, np, np)
        get_local_mass!(M, h)

        M_inv = inv(M)
        @show h
        print("M:\n")
        display(M)
        print("M_inv:\n")
        display(M_inv)
        # Stiffness matrix
        S = Array{Float64, 2}(undef, np, np)
        get_local_stiffness!(S, h)
        print("S:\n")
        display(S)
        display(-a * S)
        # Advection matrix
        A_l = M_inv * (-a .* S) + [λ 0 ; λ 0]
        print("A_l:\n")
        display(A_l)
        # Flux matrix
        F_l = Array{Float64, 2}(undef, np, 4)
        get_local_flux!(F_l, M_inv, a, α) 

        
        ## Global Matrices ##

        A_g = spzeros(Float64, np * nk, np * nk)
        get_global_advection!(A_g, A_l)
        F_g = spzeros(Float64, np * nk, np * nk)
        get_global_flux!(F_g, F_l)

        # Total Global Matrix
        G = A_g + F_g
        
        # boundary vector
        b = zeros(np * nk)
        ## Solution Vector (inital conditions) ##
        
        u_1 = zeros(t_num + 1, np * nk)
        get_inital_condition!(u_1, f, collect(x), nk)

        u_2 = copy(u_1)
        
        
        advec_local_strong_euler!(u_1, A_l, F_l, t_int, t_step, nk, g, λ)
        advec_global_strong_euler!(u_2, G, b, t_int, t_step, g, a)
        plot_DG(u_1, collect(x))
        #plot_DG(u_2, collect(x))
        #write_analytic(collect(x_an), collect(t_int), a, λ)
        #write_plot(u_2, collect(x))


        
        #=
        for t_step = .001:.001:.1
            advec_global_strong_euler!(u_1, G, b, t_int, t_step, g, a)
            #show_stability_region(Matrix(G), h, t_step)
            stab_check = stability_test(Matrix(G), h, t_step)
            @show t_step, stab_check
            print("\n\n")
            #plotDG(u_1, collect(x))
        end
        =#

    end
end
