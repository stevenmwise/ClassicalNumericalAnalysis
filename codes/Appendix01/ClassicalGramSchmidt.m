function [Q, err] = ClassicalGramSchmidt( W )
% The classical Gram-Schmidt orthogonalization process.
%
% Input
%   W(1:n,1:k) : a matrix representing a collection of k column 
%                vectors of dimension n
%
% Output
%   Q(1:n,1:k) : a collection of k orthonormal vectors of 
%                dimension n
%   err : = 0, if the columns of W are linearly independent
%         = 1, if an error has occurred
%
% WARNING: DO NOT USE THIS CODE!!
% The classical Gram-Schmidt process is numerically unstable
  n = size(W)(1);
  k = size(W)(2);
  Q = zeros(n,k);
  err = 0;
  if k > n
    err = 1;
    return;
  end
  norm_q = norm( W(:,1) );
  if norm_q > eps( norm_q )
    Q(:,1) = W(:,1)/norm_q;
  else
    err = 1;
    return;
  end
  for m=2:k
    r = W(:,m);
    for j = 1:m-1
      r = r - (W(:,m)'*Q(:,j)).*Q(:,j);
    end
    norm_q = norm( r );
    if norm_q > eps( norm_q )
       Q(:,m) = r/norm_q;
    else
       err = 1;
       return;
    end
  end
end
