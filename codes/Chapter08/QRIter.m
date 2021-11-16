function [U, EigVecs, err] = QRIter( A, maxits )
% The QR iteration method to find eigenvalues.
%
% Input
%   A(1:n,1:n) : a square matrix
%   maxits : the maximal number of iterations
%
% Output
%  U(1:n,1:n) : an approximation to an upper triangular matrix 
%               that is similar to A
%  EigVecs(1:n,1:n) : A unitary matrix whose columns are 
%                     eigenvectors of A
%  err : = 0, if no error was encountered
%        = 1, if there was an error
  err = 0;
  [U, EigVecs, err] = Hessenberg(A);
  if err == 1
    return;
  end
  for i=1:maxits
    [Q, R, err] = QRFact( U, 1 );
    if err == 1
      return;
    end
    U = R*Q;
    EigVecs = EigVecs*Q;
  end
end
