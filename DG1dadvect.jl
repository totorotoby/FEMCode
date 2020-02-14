using DelimitedFiles
using SparseArrays

function advec_global_strong_euler!(u::Array{Float64, 2}, G::SparseMatrixCSC{Float64, Int64},
                                    b::Array{Float64, 1}, t_int::StepRangeLen, step::Float64,
                                    g::Function, a::Float64)

    for (i, t) in enumerate(t_int)

        b[1] = a * g(t)
        dudt = G * u[i, :] + b

        
        u[i+1, :] = dudt .* step
        
    end

end

function get_local_mass!(M::Array{Float64, 2}, h::Float64)

    M[1,1] = 2/3 * (h/2)
    M[2,1] = 1/3 * (h/2)
    M[1, 2] = 1/3 * (h/2)
    M[2, 2] = 2/3 * (h/2)

end

function get_local_stiffness!(S::Array{Float64, 2}, h::Float64)
    
    S[1,1] = 1/2 * (h/2)
    S[2,1] = -1/2 * (h/2)
    S[1, 2] = -1/2 * (h/2)
    S[2, 2] = 1/2 * (h/2)
    
end

function get_local_flux!(F::Array{Float64, 2}, M_inv::Array{Float64, 2}, a::Float64, alpha::Float64)

    F[1,1] = M_inv[1,1]*a*alpha
    F[1,2] = -M_inv[1,1]*a*alpha 
    F[1,3] = M_inv[1,2]*a*(1 - alpha)
    F[1,4] = -M_inv[1,2]*a*(1 - alpha)
    F[2,1] = M_inv[2,1]*a*alpha
    F[2,2] = -M_inv[2,1]*a*alpha
    F[2,3] = M_inv[2,2]*a*(1-alpha)
    F[2,4] = -M_inv[2,2]*a*(1 - alpha)

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


function write_plot(u::Array{Float64, 2}, filename::String, x::Array{Float64, 1})
     
    open("plot_data", "w") do io
        writedlm(io, u)
    end
    
    x_graph = Float64[]
    for i in x
        x_graph = vcat(x_graph, [i,i])
    end
    x_graph = x_graph[2 : end - 1]
    
    open("x_data", "w") do io
        writedlm(io, x_graph)
    end
    
end


function main()

    let
        ## Parameters ##
        
        #wave-speed
        a = 1.0
        #growth-speed?
        b = 0.0
        #flux parameter
        alpha = .5
        # Left side of the interval
        l = -1
        # Right side of the interval
        r = 1
        # Number of elements
        nk = 3
        # order of approximations
        np = 2
        # Length of each element
        h = (r - l)/nk
        # elements as a range
        x = l:h:r

        # filename for plotting
        filename = "plotting_data"
        
        ## Time parameters ##
        
        # Total time
        T = 1
        # time step
        t_step = .1

        # time interval
        t_int = t_step:t_step:T
        t_num = length(t_int)
        
        # solution to graph
        u_t = Array{Float64,1}[]


        ## intial and boundary conditions ##

        f(x::Float64) = cos.(1/2*pi*x)
        g(t::Float64) = 0

        ## Local Matrices ##

        # Mass matrix
        M = Array{Float64, 2}(undef, np, np)
        get_local_mass!(M, h)

        M_inv = inv(M)

        # Stiffness matrix
        S = Array{Float64, 2}(undef, np, np)
        get_local_stiffness!(S, h)

        # Advection matrix
        A_l = a .* M_inv * S

        
        # Flux matrix
        F = Array{Float64, 2}(undef, np, 4)
        get_local_flux!(F, M_inv, a, alpha) 
        
        
        ## Global Matrices ##

        A_g = spzeros(Float64, np * nk, np * nk)
        get_global_advection!(A_g, A_l)
        F_g = spzeros(Float64, np * nk, np * nk)
        get_global_flux!(F_g, F)

        # Total Global Matrix
        G = A_g + F_g

        
        # boundary vector
        b = zeros(np * nk)
        ## Solution Vector (inital conditions) ##
        
        u = zeros(t_num + 1, np * nk)
        get_inital_condition!(u, f, collect(x), nk)

        advec_global_strong_euler!(u, G, b, t_int, t_step, g, a)

        write_plot(u, filename, collect(x))
    end
end
