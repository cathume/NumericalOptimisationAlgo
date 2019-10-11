%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%Rank one is a numerical optimisation method that is an update for
%Quasi-Newton method.
%
%----------------------------------------------------------------------------------


function [x_opt, f_opt] = rank_one(f, Df, H_init, x_init_rank1, alpha_0, alpha_1, max_iter_rank1, max_iter_secant, change_tol_secant, cahnge_tol_rank1)
%f is the objective funciton
%Df is the derivative fo the objective function
%Hk is the first guess of the approximated Hessian matrix
%alpha_0 and alpha_1 are first guesses of the step size
%max_iter... is the max iterations
%change_tol... is the tolerance for difference between iterations
iter = 1;
x_curr = x_init_rank1;
Hk = H_init;
while iter <= max_iter_rank1
    if Df(x_curr)==0
        fprintf('You have the optimum solution %f', x_curr);
        break
    end
    d = -Hk*Df(x_curr);
    Dphi =@(alpha) Df(x_curr+alpha*d)'*d;
    [alpha] = secant_line_search(alpha_0,alpha_1,Dphi,max_iter_secant,change_tol_secant);
    x = x_curr+ alpha*d;
    delta_x = x-x_curr;
    delta_g= Df(x)-Df(x_curr);
    Hk = Hk + ((delta_x - Hk*delta_g)*(delta_x-Hk*delta_g)')/(delta_g'*(delta_x-Hk*delta_g));
    x_curr = x;
    disp('iteration number');
    disp(iter);
    disp('x');
    disp(x_curr);
    norm_2 = norm([1;1]-x_curr);
    disp('norm');
    disp(norm_2);
    disp('g');
    disp(norm(Df(x_curr)));
    disp('f');
    disp(f(x_curr));
    iter = iter+1;
    
end        
x_opt = x_curr;
f_opt = f(x_opt);