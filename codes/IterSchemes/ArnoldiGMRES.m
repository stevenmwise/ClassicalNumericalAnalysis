function [Q, H, err] = ArnoldiGMRES( A, r, m )
% The Arnoldi algorithm to compute an orthonormal basis of
%   K_k(A, r) = span{ r, Ar, ..., A^{m-1} r }
% the Krylov subspace of order m. This will be eventually used in GMRES
%
% Input:
%   A(1:n,1:n) : the system matrix
%   r(1:n) : the initial residual
%   m : the size of the Krylov subspace
%
% Output:
%  Q(1:n, 1:m) : the columns of this matrix are the orthonormal basis
%  H(1:m, 1:m-1) : the upper Hessenberg matrix that is defined by
%    AQ_k = Q_{k+1} H
%  with Q_i being the first columns of Q
%  err : = 0, if the algorithm proceeded to completion 
%        = 1, if H(j+1,j) = 0 at some point
  err = 0;
  n = size(A,1);
  Q = zeros(n,m);
  H = zeros(m,m-1);
  norm_r = norm( r );
  if norm_r < eps( norm_r )
    err = 1;
    return;
  end
  Q(:,1) = r/norm_r;
  for j=1:m-1
    Q(:,j+1) = A*Q(:,j);
    for i=1:j
      H(i,j) = Q(:,i)'*Q(:,j+1);
      Q(:,j+1) = Q(:,j+1) - H(i,j)*Q(:,i);
    end
    norm_r = norm( Q(:,j+1) );
    if norm_r < eps( norm_r )
      err = 1;
      return;
    end
    H(j+1,j) = norm_r;
    Q(:,j+1) = Q(:,j+1)/norm_r;
  end
end
