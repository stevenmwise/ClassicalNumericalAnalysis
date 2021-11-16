function [root,count,err] = NewtonRoot(xin, p, q, tol, maxits)
%
% This function calculates the pth root of q ,
%
%   q^(1/p)
%
% using Newton's method. For simplicity, we assume that p is a 
% positive integer and q is a positive real number
%
% Input:
%   xin : initial guess
%   p : the positive integer degree of the root
%   q : the positive number whose pth root is to estimated
%   tol : stopping tolerance
%   maxits : the maximal number of iterations
%
% Output:
%   root : approximation of the root
%   count : number of newton iterations required to compute the 
%           root
%   err: = 0, if the algorithm proceeded to completion
%        = 1, if an error was encountered
%
  root = NaN;
  count = 0;
  err = 0;
  if int32(p) ~= p || p < 0
    disp('Error: p must be a positive integer'); 
    err = 1;
  return
  end
  if q < 0
    disp('Error: q must be positive');
    err = 1;
  return
  end
  diff = 1.0;
  x = xin;
  while diff > tol && count < maxits
    xo = x;
    x = xo - fn(xo,p,q)/dfn(xo,p,q);
    diff = abs(x-xo);
    count = count+1;
  end
  if count >= maxits
    err = 1;
  end
  root = x;
end

function y = fn(x,p,q)
  y = x^p-q;
end

function y = dfn(x,p,q)
  y = p*x^(p-1);
end
