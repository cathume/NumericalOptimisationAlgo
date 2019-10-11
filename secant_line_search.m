%----------------------------------------------------------------------------------
%Author: Catherine Hume
%
%Secant line search is used in one dimensional optimisation. Usually, I use
%this to find the optimal step size in numerical optimisation methods.
%----------------------------------------------------------------------------------


function [alpha] = secant_line_search(alpha_0, alpha_1, Dphi, max_iter_secant, change_tol_secant)
%alpha_0 and alpha_1 are two inital values
%Dphi is is the derivative of the objective function
%max_iter is the maximum number of iterations
%change_tol is the change tolerance
iter = 1;
while iter<= max_iter_secant
   alpha = (alpha_0*Dphi(alpha_1) - alpha_1*Dphi(alpha_0))/(Dphi(alpha_1)-Dphi(alpha_0));
  if abs(Dphi(alpha))< change_tol_secant
      %or abs((alpha-alpha_1)/alpha) < change_tol
      %make sure to use absolute value!
       break
  end
   alpha_0 = alpha_1;
   alpha_1= alpha;
   iter = iter+1;
end
%display your output by uncommenting the line below. 
% disp(alpha)
