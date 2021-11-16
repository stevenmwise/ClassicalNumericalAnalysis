% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.  
%          If v is complex, delete "real" commands.

  function w = chebfft(v)
  N = length(v)-1; if N==0, w=0; return, end
  x = cos((0:N)'*pi/N);
  ii = 0:N-1;
  v = v(:); V = [v; flipud(v(2:N))];      % transform x -> theta          
  U = real(fft(V));
  W = real(ifft(1i*[ii 0 1-N:-1]'.*U));
  w = zeros(N+1,1);
  w(2:N) = -W(2:N)./sqrt(1-x(2:N).^2);    % transform theta -> x     
  w(1) = sum(ii'.^2.*U(ii+1))/N + .5*N*U(N+1);     
  w(N+1) = sum((-1).^(ii+1)'.*ii'.^2.*U(ii+1))/N + ...
              .5*(-1)^(N+1)*N*U(N+1);
