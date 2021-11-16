function [R, err] = CholeskyDecomposition( A )
% Choleksy factorization of the HPD matrix A:
%   A = R'*R
%
% Input:
%   A(1:n,1:n) : an HPD matrix
%
% Output:
%   R(1:n,1:n) : The upper triangular matrix satisfying A = R'*R
%   err : = 0 if the decomposition finished successfully
%         = 1 if an error occurred
  err = 0;
  n = size(A,1);
  R = zeros(n,n);
  R(1,1) = sqrt( A(1,1) );
  for i=2:n
    for j=1:i-1
      sum = A(i,j);
      for k=1:j-1
        sum = sum - R(k,i)* conj( R(k,j) );
      end
      if abs( R(j,j) ) > eps( R(j,j) )
        R(j,i) = sum/R(j,j);
      else
        err = 1;
        return;
      end
    end
    sum = A(i,i);
    for k=1:i-1
      sum = sum - R(k,i)*conj( R(k,i) );
    end
    R(i,i) = sqrt( sum );
  end
end
