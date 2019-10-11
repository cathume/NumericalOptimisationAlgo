%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%FISTA refers to the newer "fast iterative shinkrage thresholding
%algorithm" (2009) Beck and Teboulle.
%coded here independently for a research project
%----------------------------------------------------------------------------------


function [x_opt,f_opt, f_val,cpu_time] = FISTA (x_init, max_iter, f, Df, L, change_tol)
%x_init is the initial guess for an x-value
% max_iter is the maximum number of iterations
% f is the objective function
%Df is the derivative of the objective function
%L is a Lipschitz constant of A
%f_val is a vector storing the output values of the function
x_curr = x_init;
y = x_init;
t_curr = 1;
iter = 1;
f_val = zeros(max_iter,1);
cpu_time = zeros(max_iter,1);
tic;
while iter < max_iter
    
    cpu_time (iter) = toc;
    x = y-(Df(y).*(1/L));
    t = (1 + sqrt(1 + 4* (t_curr)^2))/2;
    y = x + ((t_curr -1)/t)*(x-x_curr);
    if (norm(x_curr - x))< change_tol
        break;
    end
    t_curr = t;
    x_curr = x;
    disp("iter = " + iter);
    iter = iter + 1;
    disp("x_curr = ");
    disp(x_curr);
    f_val(iter) = f(x_curr);
end
x_opt = x_curr;
disp("x_opt = ");
disp(x_opt);
f_opt = f(x_opt);
disp("f_opt = ");
disp(f_opt);
return