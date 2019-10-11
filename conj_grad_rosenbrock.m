function [x_opt,f_opt] = conj_grad_rosenbrock(f, Df, Hf, x_init, max_iter, change_tol_conj_grad_rosenbrock)
%this prgram optimises the objective function with the conjugate gradient
%when we have a function that isnt of standard quadratic form! The Hessian
%matrix is calculated at every iteration. Feel free to shoot yourself
%before coding anymore :p



%gf is the gradient of the objective function
%Hf is the hessian matric of the objective function
%x_init is the first guess
%max_iter is maximum number of iterations
%relchng_tol is the relative change tolerance
iter = 1;
x_curr = x_init;
g = Df(x_curr);
h = Hf(x_curr);
if g == 0
    disp('Your first guess is the optimal value!');
else
    d = -g;
    while iter<= max_iter
        alpha = -1*((g'*d)/(d'*h*d));
        x = x_curr + alpha*d;
        g = Df(x);
        h = Hf(x);
        if g==0
            break
        end
        beta = (g'*h*d)/(d'*h*d);
        d = beta*d-g;
        if norm(x-x_curr)/norm(x_curr) < change_tol_conj_grad_rosenbrock
            break
        end
        x_curr = x;
        disp('iteration number:');
        disp(iter);
        disp(x_curr);
        disp(f(x_curr));
        iter = iter+1;
    end
    x_opt = x;
    f_opt = f(x_opt);
    
end
disp(x_opt);
disp(f_opt);