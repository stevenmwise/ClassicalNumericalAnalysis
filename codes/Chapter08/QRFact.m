function [Q, R, err] = QRFact( A, computeQ )
% The QR factorization using Householder reflectors.
%
% Input
%   A(1:m,1:n) : a rectangular matrix representing a collection 
%                of n column vectors of dimension m
%   computeQ : = 1, the matrix Q is computed
%              = 0, the matrix Q is not computed. The algorithm 
%                returns the identity
%
% Output
%   Q(1:m,1:m) : a unitary matrix of dimension n
%   R(1:m,1:n) : an rectangular upper triangular square matrix 
%                such that A = QR
%   err : = 0, if the columns of A are linearly independent
%         = 1, if an error has occurred
  m = size(A,1);
  n = size(A,2);
  Q = eye(m,m);
  R = A;
  err = 0;
  if n>m
    err = 1;
    return;
  end
  for k=1:n
    nn = norm( R(k:m,k) );
    e = zeros(m-k+1,1);
    e(1) = 1;
    v = nn*e + R(k:m,k);
    norm_v = norm( v );
    if norm_v > eps( norm_v )
      w = (1.0/norm_v)*v;
    else
      err = 1;
      return;
    end
    R(k:m,k:n) = R(k:m,k:n) - 2.0*(w*w')*R(k:m,k:n);
    if computeQ > 0
      for t=1:m
        Q(k:m,1:m) = Q(k:m,1:m) - 2.0*(w*w')*Q(k:m,1:m);
      end
    end    
  end
  if computeQ > 0
    Q = Q';
  end
end
