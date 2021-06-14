function [Q, R, err] = ModifiedGramSchmidt( A )
% The modified Gram-Schmidt orthogonalization process.
%
% Input
%   A(1:m,1:n) : a matrix representing a collection of n column vectors of dimension m
%
% Output
%   Q(1:m,1:n) : a collection of k orthonormal vectors of dimension n
%   R(1:n,1:n) : an upper triangular square matrix such that A = QR
%   err : = 0, if the columns of A are linearly independent
%         = 1, if an error has occurred
  m = size(A)(1);
  n = size(A)(2);
  Q = zeros(m,n);
  R = zeros(n,n);
  err = 0;
  if n>m
    err = 1;
    return;
  end
  V = A;
  for i=1:n
    R(i,i) = norm( V(:,i) );
    if R(i,i) > eps( R(i,i) );
      Q(:,i) = (1./R(i,i)) * V(:,i);
    else 
      err = 1;
      return;
    end
    for j = i+1:n
      R(i,j) = V(:,j)'*Q(:,i);
      V(:,j) = V(:,j) - R(i,j)*Q(:,i);
    end
  end
end
