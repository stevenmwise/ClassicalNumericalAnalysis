function [x, its, err] = MinRes( A, x0, f, maxit, tol )
% The method of minimal residuals to approximate the solution to
%
%   Ax = f
%
% with A HPD.
%
% Input
%   A(1:n,1:n) : the system matrix,
%   x0(1:n) : the initial guess 
%   f(1:n) : the right hand side vector
%   maxit : the maximal number of iterations
%   tol : the tolerance
%
% Output
%   x(1:n) : the approximate solution to the linear system of equations
%   its : the number of iterations
%   err : = 0, if the tolerance is reached in less than maxit iterations
%         = 1, if the tolerance is not reached
  err = 0;
  x = x0;
 
  for its=1:maxit
     r = f-A*x;
     p = A*r;
     p_norm = norm( p );
     if p_norm < tol 
       return;
     end
     alpha = r'*p/(p_norm*p_norm);
     x = x + alpha*r;
  end
  err = 1;
end
