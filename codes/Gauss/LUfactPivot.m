function [L, U, P, err] = LUfactPivot( A )
% LU factorization with pivoting of a square matrix.
%
% Input
%   A(1:n,1:n) : the matrix to be factorized
%
% Output
%   L(1:n,1:n), U(1:n,1:n), P(1:n,1:n) : the factors in A = P'LU
%   err : = 0 if no error was encountered
%         = 1 if a division by zero occurred
  n = size(A,1);
  U = A;
  L = eye(n);
  P = eye(n);
  err = 0;
  for k=1:n-1
    i=k;
    for t=k:n
      if abs( U(t,k) ) > abs( U(i,k) )
        i=t;
      end
    end
    for t=k:n
      temp = U(k,t);
      U(k,t) = U(i,t);
      U(i,t) = temp;
    end
    for t=1:k-1
      temp = L(k,t);
      L(k,t) = L(i,t);
      L(i,t) = temp;
    end
    for t=1:n
      temp = P(k,t);
      P(k,t) = P(i,t);
      P(i,t) = temp;
    end
    if abs( U(k,k) ) > eps( U(k,k) )
      for j=k+1:n
        L(j,k) = U(j,k)/U(k,k);
        for t=k:n
          U(j,t) = U(j,t) - L(j,k)*U(k,t);
        end
      end
    else
      err = 1;
      return;
    end
  end
end
  
