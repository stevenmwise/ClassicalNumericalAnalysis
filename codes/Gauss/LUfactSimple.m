function [L, U, err] = LUfactSimple( A )
% LU factorization of a square matrix.
%
% Input
%   A(1:n,1:n) : the matrix to be factorized
%
% Output
%   L(1:n,1:n), U(1:n,1:n) : the factors in A = LU, if Gaussian elimination proceeds to completion
%   err : = 0 if Gaussian elimination proceeds to completion
%         = 1 if a zero pivot is encountered
  n = size(A,1);
  U = A;
  L = eye(n);
  err = 0;
  for k=1:n-1
    for j=k+1:n
      if abs( U(k,k) ) > eps( U(k,k) )
        L(j,k) = U(j,k)/U(k,k);
      else
        err = 1;
        return;
      end
      for t = k:n
        U(j,t) = U(j,t)-L(j,k)*U(k,t);
      end
    end
  end
end
