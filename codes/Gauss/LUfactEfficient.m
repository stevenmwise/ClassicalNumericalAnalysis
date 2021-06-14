function [Afact, swaps, err] = LUfactEfficient( A )
% Efficient implementation of LU factorization with partial pivoting.
%
% Input
%   A(1:n,1:n) : the matrix to be factorized into A = P'LU
%
% Output
%   Afact(1:n,1:n) : a matrix containing the matrices L and U below and above the diagonal, respectively
%   swaps : the indices that indicate the permutations (row swaps)
%   err : = 0 if no the factorization finished successfully
%         = 1 if a division by zero was encountered
  n=size(A,1);
  swaps = 1:n;
  err = 0;
  for k=1:n-1
    i = k;
    for t=k:n
      if abs( A(t,k) ) > abs( A(i,k) )
        i=t;
      end
    end
    tt = swaps(k);
    swaps(k) = swaps(i);
    swaps(i) = tt;
    if abs( A( swaps(k), k ) <= eps( A( swaps(k), k) )
      err = 1;
      return;
    end
    for i=k+1:n 
      xmult = A(swaps(i),k)/A(swaps(k),k);
      A( swaps(i), k ) = xmult;
      for j=k+1:n
        A(swaps(i),j) = A(swaps(i),j) - xmult*A(swaps(k),j);
      end
    end
  end
  Afact = A;
end
