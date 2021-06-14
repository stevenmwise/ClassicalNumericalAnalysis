function [U, EigVecs, err] = QRshifts( A, maxits )
% The QR iteration method with shifts to find eigenvalues.
%
% Input
%   A(1:n,1:n) : a square matrix
%   maxits : the maximal number of iterations
%
% Output
%  U(1:n,1:n) : an approximation to an upper triangular matrix that is similar to A
%  err : = 0, if no error was encountered
%        = 1, if there was an error
  err = 0;
  [U, EigVecs, err] = Hessenberg(A);
  if err == 1
    return;
  end
  n = size(A,1);
  id = eye(n,n);
  for i=1:maxits
    mu = U(n,n);
    [Q, R, err] = QRtriangularization( U - mu*id, 1 );
    if err == 1
      return;
    end
    U = R*Q + mu*id;
    
    EigVecs = EigVecs*Q;
  end
end
