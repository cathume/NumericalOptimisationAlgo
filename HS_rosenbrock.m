function [x_opt,f_opt,f_val_CG] = HS_rosenbrock(f, Df, x_init, max_iter, change_tol_conj_grad_rosenbrock, alpha_0, alpha_1, max_iter_secant, change_tol_secant)
%gf is the gradient of the objective function
%Hf is the hessian matric of the objective function
%x_init is the first guess
%max_iter is maximum number of iterations
%relchng_tol is the relative change tolerance
iter = 1;
x_curr = x_init;
g = Df(x_curr);
if g == 0
    disp('Your first guess is the optimal value!');
else
    d = -g;
    while iter<= max_iter
        Dphi =@(alpha) Df(x_curr+alpha*d)'*d;
        [alpha] = secant_line_search(alpha_0,alpha_1,Dphi,max_iter_secant,change_tol_secant);
        disp(alpha);
        x = x_curr + alpha*d;
        g = Df(x);
        if g==0
            break
        end
        beta = (Df(x)'*(Df(x)-Df(x_curr))/(d'*(Df(x)-Df(x_curr))));
        d = beta*d-g;
        if norm(x-x_curr)/norm(x_curr) < change_tol_conj_grad_rosenbrock
            break
        end
        x_curr = x;
        disp('iteration number:');
        disp(iter);
        disp('x_curr:');
        disp(x_curr);
        disp('g');
        disp(norm(g));
        disp ('f(x)');
        disp(f(x_curr));
        iter = iter+1;
    end
    x_opt = x;
    f_opt = f(x_opt);
    
end
%disp(x_opt);
%disp(f_opt);