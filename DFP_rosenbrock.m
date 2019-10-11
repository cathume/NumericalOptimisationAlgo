%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%DFP refers to the Davidson, Fletcher, and Powell update to the
%approximation of the inverse Hessian matrix from Newton's method. Unlike
%rank-one it assures that the approximation will be positive definite
%
%----------------------------------------------------------------------------------



function [x_opt, f_opt,f_val_DFP] = DFP_rosenbrock(f, Df, H_init, x_init_DFP, alpha_0, alpha_1, max_iter_DFP, max_iter_secant, change_tol_secant, change_tol_DFP, xstar)
%f is the objective funciton
%Df is the derivative fo the objective function
%Hk is the first guess of the approximated Hessian matrix
%alpha_0 and alpha_1 are first guess of the step size
%max_iter... is the max iterations
%change_tol... is the tolerance for difference between iterations
iter = 1;
x_curr = x_init_DFP;
Hk = H_init;
f_val_DFP = zeros(max_iter_DFP,1);
normie_DFP = zeros(max_iter_DFP,1);
while iter <= max_iter_DFP
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
    Hk = Hk + ((delta_x*delta_x')/(delta_x'*delta_g))-(((Hk*delta_g)*(Hk*delta_g)')/(delta_g'*Hk*delta_g));
    if norm(x-x_curr/x)< change_tol_DFP
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
    iter = iter+1;
    f_val_DFP(iter) = f(x_curr);
    normie_DFP(iter) = norm(x_curr-xstar);
    
end        
x_opt = x_curr;
f_opt = f(x_opt);