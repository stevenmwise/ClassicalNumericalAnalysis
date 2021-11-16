function [x, err] = Solve( A, f )
% Solves the system Ax = f.
%
% Input
%   A(1:n,1:n) : the coefficient matrix
%   f(1:n) : the RHS vector
%
% Output:
%   x(1:n) : the solution vector
%   err : = 0 if the solution was found successfully
%         = 1 if an error occurred during the process
  n = length(f);
  err = 0;
  [AA, swaps, err] = LUFactEfficient( A );
  if err == 1
    return
  end
  for k=1:n-1
    for i=k+1:n
      f( swaps(i) ) = f( swaps(i) ) - AA( swaps(i), k ) ...
        *f( swaps(k) );
    end
  end
  if abs( AA( swaps(n), n ) ) > eps( AA(swaps(n),n) )
    x(n) = f( swaps(n) )/AA( swaps(n), n );
  else
    err = 1;
    return;
  end
  for i=n-1:-1:1
    xsum = f( swaps(i) );
    for j=i+1:n
      xsum = xsum - AA( swaps(i), j )*x(j);
    end
    if abs( AA( swaps(i), i ) ) > eps( AA( swaps(i), i ) )
      x(i) = xsum/AA( swaps(i), i );
    else
      err = 1;
      return;
    end
  end
end
