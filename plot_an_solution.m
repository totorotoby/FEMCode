

an_sol = animatedline;
num_sol = animatedline;
axis([-1, 1. -.5 1]);
an_data = readmatrix("an_data");
x_an_axis = readmatrix("x_an_axis");
x_num_axis = readmatrix("x_num_data");
num_data = readmatrix("num_data");

display(num_data);


for t_step = 1:size(an_data,1)
    
    clearpoints(an_sol)
    addpoints(an_sol, x_an_axis, an_data(t_step, :));
    pause(.1);
end