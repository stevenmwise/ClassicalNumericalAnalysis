function [x, err] = TriDiagonal( a, b, c, f )
% Solution of a linear system of equations with a tridiagonal 
% matrix.
%
%   a(k) x(k-1) + b(k) x(k) + c(k) x(k+1) = f(k)
%
% with a(1) = c(n) = 0.
%
% Input
%   a(1:n), b(1:n), c(1:n) : the coefficients of the system 
%                            matrix
%   f(1:n) : the right hand side vector
%
% Output
%   x(1:n) : the solution to the linear system of equations, if 
%            no division by zero occurs
%   err : = 0, if no division by zero occurs
%         = 1, if division by zero is encountered
  n = length(f);
  alpha = zeros(n,1);
  beta = zeros(n,1);
  err = 0;
  if abs(b(1)) > eps( b(1) )
    alpha(1) = -c(1)/b(1);
    beta(1) = f(1)/b(1);
  else
    err = 1;
    return;
  end
  for k=2:n
    denominator = a(k)*alpha(k-1) + b(k);
    if abs(denominator) > eps( denominator )
      alpha(k) = - c(k)/denominator;
      beta(k) = ( f(k) - a(k)*beta(k-1) )/denominator;
    else
      err = 1;
      return;
    end
  end
  if abs(a(n)*alpha(n-1) + b(n)) > eps( b(n) )
     x(n) = ( f(n) - a(n)*beta(n-1) )/( a(n)*alpha(n-1) + b(n) );
  else
    err = 1;
    return;
  end
  for k=n-1:-1:1
    x(k) = alpha(k)*x(k+1) + beta(k);
  end
end
