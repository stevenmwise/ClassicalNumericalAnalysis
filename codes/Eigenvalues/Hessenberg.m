function [M, Q, err] = Hessenberg(A)
% The following algorithm reduces the matrix A to Hessenberg form
%
% M = Q'AQ
%
% Input
%   A(1:n,1:n) : a square matrix
%
% Output
%  M(1:n,1:n) : a lower Hessenberg matrix that is similar to A
%  Q(1:n,1:n) : the unitary matrix that realizes the transformation
%  err : = 0, if no error was encountered
%        = 1, if there was an error
  [m,n] = size(A);
  M = A;
  err = 0;
  if m ~= n
    err = 1;
    return;
  end
  Q = eye(n);
  for k=1:n-2
    nn = norm( M(k+1:n,k) );
    e = zeros(n-k,1);
    e(1) = 1;
    v = nn*e + M(k+1:n,k);
    norm_v = norm(v);
    if norm_v < eps( norm_v )
      err = 1;
      return;
    end
    v = (1.0/norm_v)*v;
    M(k+1:n,k:n) = M(k+1:n,k:n) - 2.0*(v*v')*M(k+1:n,k:n);
    M(1:n,k+1:n) = M(1:n,k+1:n) - 2.0*M(1:n,k+1:n)*(v*v');
    Q(1:n,k+1:n) = Q(1:n,k+1:n) - 2.0*Q(1:n,k+1:n)*(v*v');
  end
end
