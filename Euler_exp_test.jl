using Plots

l = 0

r = 0

b = 4

t_step = .01
time = t_step : t_step : 1
x_plot = 0 : t_step : 1
us = zeros(length(x_plot))
us[1] = 1
for (i, t) in enumerate(time)

    us[i+1] = us[i] + t_step * ( b * us[i])
    
end

plt = Plots.plot()
plot!(x_plot, us)
