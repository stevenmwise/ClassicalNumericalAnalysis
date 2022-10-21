function [x, err] = CyclicallyTriDiagonal( a, b, c, f )
% Solution of a cyclically tridiagonal linear system of equations:
%
%   b(1) x(  1) + c(1) x(  2) + a(1) x(  n) = f(1),
%   a(k) x(k-1) + b(k) x(  k) + c(k) x(k+1) = f(k),
%   c(n) x(  1) + a(n) x(n-1) + b(n) x(  n) = f(n).
%
% Input
%   a(1:n), b(1:n), c(1:n) : the coefficients of the system 
%                            matrix
%   f(1:n) : the right hand side vector
%
% Output
%   x(1:n) : the solution to the linear system of equations, if
%            no division by zero is encountered
%   err : = 0, if division by zero does not occur
%         = 1, if division by zero occurs
  n = length(f);
  err = 0;
  newa = a(2:n);
  newa(1) = 0;
  newb = b(2:n);
  newc = c(2:n);
  newc(n-1) = 0;
  newf = f(2:n);
  [u, err] = TriDiagonal( newa, newb, newc, newf );
  if err == 1
    return;
  end
  newf = zeros(n-1, 1);
  newf(1) = - a(2);
  newf(n-1) = -c(n);
  [v, err] = TriDiagonal( newa, newb, newc, newf );
  if err == 1
    return;
  end
  x = zeros(n,1);
  denominator = b(1) + c(1)*v(1) + a(1)*v(n-1);
  if abs(denominator) > eps( denominator )
    x(1) =  ( f(1) - c(1)*u(1) - a(1)*u(n-1) )/denominator;
  else
    err = 1;
    return;
  end
  for k=2:n
    x(k) = u(k-1) + x(1)*v(k-1);
  end 
end
