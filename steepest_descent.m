%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%Steepest descent optimisation method.
%
%----------------------------------------------------------------------------------


function [x_opt] = steepest_descent (Df, x_init, max_iter, change_tol,alpha_0,alpha_1,Dphi, max_iter_secant,change_tol_secant);
%Df is the gradient of the objective function
%x_init is the initla x value you want to test
%max_iter is the maximum number of iterations
%change_tol is the change tolerance
%all elements from alpha_0 and onwards are defined in secant method code

iter = 1;
x_curr = x_init;
while iter<= max_iter
    g = Df(x_curr);
    if g == 0
%         disp('This is your optimised value!');
%         disp (x_opt);
        break
    end
    [alpha] = secant(alpha_0, alpha_1, Dphi, max_iter_secant, change_tol_secant);
    x = x_curr - g* alpha;
    if abs((x-x_curr)/x)< change_tol
        break
    end
    x_curr = x;
    iter = iter + 1;

end

if iter == max_iter
    disp('Max number of iterations reached.');
end
x_opt = x;
disp(x_opt);