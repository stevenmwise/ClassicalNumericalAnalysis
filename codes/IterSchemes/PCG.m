function [x, its, err] = PCG( A, x0, f, Binv, maxit, tol )
% The preconditioned conjugate gradient method to approximate the solution to
%
%   Ax = f
%
% with A HPD.
%
% Input
%   A(1:n,1:n) : the system matrix
%   x0(1:n) : the initial guess 
%   f(1:n) : the right hand side vector
%   Binv(1:n,1:n) : the inverse of the preconditioner
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
  r = f - A*x;
  p = Binv*r;
  z = p;
  initerror = sqrt( r'*p );
  for its=1:maxit
     Ap = A*p;
     denom = p'*Ap;
     if denom < tol
       err = 0;
       return;
     end
     theta = (1./denom)*(z'*p);
     x = x + theta*p;
     r = r - theta*Ap;
     z = Binv*r;
     nu = (1./denom)*(z'*p);
     p = z + nu*p;
     if sqrt( r'*z )/initerror < tol 
       err = 0;
       return;
     end
  end
  err = 1;
end
