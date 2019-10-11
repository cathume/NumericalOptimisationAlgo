%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%basic algorithm for conjugate gradient method
%----------------------------------------------------------------------------------


function [x_opt_2, f_val_CG,normie_CG] = conj_grad(f, Df, Hk, x_init, max_iter, change_tol_conj_grad, xstar)
%gf is the gradient of the objective function
%Hf is the hessian matric of the objective function
%x_init is the first guess
%max_iter is maximum number of iterations
%relchng_tol is the relative change tolerance
iter = 1;
x_curr = x_init;
g = Df(x_curr);
f_val_CG = zeros(max_iter, 1);
normie_CG = zeros(max_iter,1);
if g == 0
    disp('Your first guess is the optimal value!');
else
    d = -g;
    while iter <= max_iter 
        g = Df(x_curr);
        alpha = -1 * (((g)'*d) / ((d)'* Hk * d));
        x = x_curr + alpha*d;
        g = Df(x_curr);
        if g==0
            break
        end
        beta = (g'*Hk*d)/(d'*Hk*d);
        d = beta*d-g;
        if norm(x-x_curr)/norm(x_curr) < change_tol_conj_grad
            break
        end
        x_curr = x;
        disp('iter');
        disp(iter);
        iter = iter+1;
        disp('x');
        disp(x_curr);
        disp('beta');
        disp(beta);
        disp('direction');
        disp(d);
        f_val_CG(iter) = f(x_curr);
        normie_CG(iter) = norm(x_curr - xstar);
    end
    x_opt_2 = x;
  
    
end
%disp(x_opt);