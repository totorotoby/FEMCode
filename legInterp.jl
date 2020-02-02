#=



=#
using Plots

#=

function to plot interpolatent of arbitrary order on [-1, 1] with Legendre Polynomials

=#

function plotU_h(u_coords::Array{Float64}, plot_points::Array{Float64})
    plot()
    y_plot = u_coord[1].*phi_1(plot_points) .+  u_coord[2].*phi_2(plot_points) #.+  u_coord[2].*phi_3(plot_points)
    plot!(plot_points, u(plot_points))
    return plot!(plot_points, y_plot)
end

# number of Legendre Polynomials
O = 2

# number of interior interpolation points in a single element
in_step = 2/(O-1)

# random u to test things out with
u(x) = x.^5 .+3x .- 1

phi_1(x) = 1/sqrt(2)
phi_2(x) = sqrt(3/2).*x
phi_3(x) = 3*sqrt(5)/2*sqrt(2).*x.^2 .- sqrt(5)/2*sqrt(2)


phi = [phi_1, phi_2]

# interpolation points on the unit interval [-1,1]
xi = -1:in_step:1

# vandermonde matrix
V = zeros(O,O)

for i in 1:O
    for j in 1:O
        V[i,j] = phi[j](xi[i])
    end
end

# points on u that we are trying to find linear combinations of basis functions to equal.
u_v = [u(xi[j]) for j in 1:O]

# coordinates for basis functions
u_coord = V\u_v

#### plotting u and combinations of phi ####
plot_step = .05
plot_points = -1:plot_step:1
#@show(typeof(u_coord), typeof(collect(plot_points)))
plotU_h(u_coord, collect(plot_points))


