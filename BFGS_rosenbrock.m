%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%BFGS (Broyden, Flethcer, Goldfarb, Shannon) is an update for the
%quasi-newton method approximation of the Hessian matrif of the objective
%function.
%----------------------------------------------------------------------------------

function [x_opt_1, f_opt,f_val_BFGS] = BFGS_rosenbrock(f, Df, H_init, x_init_BFGS, alpha_0, alpha_1, max_iter_BFGS, max_iter_secant, change_tol_secant, change_tol_BFGS)
%f is the objective funciton
%Df is the derivative fo the objective function
%Hk is the first guess of the approximated Hessian matrix
%alpha_0 and alpha_1 are first guess of the step size
%max_iter... is the max iterations
%change_tol... is the tolerance for difference between iterations
iter = 1;
x_curr = x_init_BFGS;
Hk = H_init;
f_val_BFGS = zeros(max_iter_BFGS, 1);
while iter <= max_iter_BFGS
    if Df(x_curr)==0
        fprintf('You have the optimum solution %f', x_curr);
        break
    end
    d = -1*Hk*Df(x_curr);
    Dphi =@(alpha) Df(x_curr+alpha*d)'*d;
    [alpha] = secant_line_search(alpha_0,alpha_1,Dphi,max_iter_secant,change_tol_secant);
    x = x_curr+ alpha*d;
    delta_x = (x-x_curr);
    delta_g= (Df(x)-Df(x_curr));
    Hk = Hk + (((delta_g)*(delta_g)')/((delta_g)'*(delta_x))) - ((Hk*(delta_x)*(delta_x)'*Hk)/((delta_x)'*Hk*(delta_x)));
    if norm(x-x_curr/x)< change_tol_BFGS
        break
    end
    x_curr = x;
    disp('iteration number');
    disp(iter);
    disp('x');
    disp(x_curr);
    disp('g');
    disp(norm(Df(x_curr)));
    disp('f');
    disp(f(x_curr));
    f_val_BFGS(iter) = f(x_curr);
    iter = iter+1;
    
end        
x_opt_1 = x_curr;
f_opt = f(x_opt_1);