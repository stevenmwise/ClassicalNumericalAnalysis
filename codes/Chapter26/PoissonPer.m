function [ePSpec,eFDiff] = PoissonPer(L, N)
%
% This function computes pseudo-spectral and finite difference 
% approximations to the periodic Poisson equation:
%
%   -D_{xx} u = f, where u and f are [0,L]-periodic
%
% using the Fast Fourier Transform (FFT) algorithm.
%
% Input
%
%   L : the size of the domain
%   N : the number of points in the periodic grid. Generally 
%       this should be a positive, even number.
%
% Output
%
%   ePSpec : the max norm of the error in the pseudo-spectral 
%          computation.
%   eFDiff : the max norm of the error in the finite difference 
%          computation.
%
% h is the grid spacing; x holds the nodal values of the grid:
  h = L/N;
  x = (1:N)*h;
%
% Fine grid parameters for representing the exact solution:
  Nf = 256;
  hf = L/Nf;
  xf = (1:Nf)*hf;
%
% q is a useful factor that allows us to easily change the 
%   domain size:
  q = 2*pi/L;
%
% Wave numbers, k: This definition is really convenient,
% essentially making k periodic. Note that our array index
% starts at 1, as usual in MATLAB, but the value of k starts
% at 0:
  k = [0:N/2-1 N/2 -N/2+1:-1];
%
% Right-Hand-Side function, f, forced to be discrete mean-zero:
  f = q*q*(cos(q*x).*cos(q*x)-sin(q*x)).*exp(sin(q*x));
  f = -(f-h*sum(f)/L);
%
% Exact solution, forced to be discrete mean-zero:
  uExact = exp(sin(q*x));
  uExactMass = h*sum(uExact)/L;
  uExact = uExact-uExactMass;
%
% Eigenvalues of the pseudo-spectral derivative operator -D^2:
  eigenPSpec = q*q*k.*k;
%
% Modified to avoid dividing by zero below:
  eigenPSpec(1) = 1.0;
%
% Eigenvalues of the finite difference operator -\Delta_h. 
% Note that we have applied a shift here to be consistent with 
% the pseudo-spectral case:
  eigenFDiff = 4.0*sin(q*(x-h)/2).*sin(q*(x-h)/2)/(h*h);
%
% Modified to avoid dividing by zero below:
  eigenFDiff(1) = 1.0;
%
% Approximations: Note that MATLAB fft's do not respect the fact
% that our data points are real numbers. So we need to get rid 
% of any imaginary numbers in the results:
  uPSpec = real((ifft(fft(f)./eigenPSpec)));
  uFDiff = real((ifft(fft(f)./eigenFDiff)));
%
% Error computations:
  ePSpec = max(abs(uPSpec-uExact));
  eFDiff = max(abs(uFDiff-uExact));
%
% Redefine uExact for plotting:
  uExact = exp(sin(q*xf));
  uExact = uExact-uExactMass;
  
  hf = figure(1)

  subplot(2,1,1);
  plot(x,uPSpec,'s',x,uFDiff,'o',xf,uExact,'k-')
  grid on, xlabel x, ylabel 'exact and approximate solutions';
  title(['Pseudo-Spectral and Finite Difference ', ...
    'Approximations: N = ', num2str(N)]);
  axis([0,L,-1.1,1.7])
  set(gca,'xTick',0:L/8:L)
  text(1.25,0.425,['max error ps = ' ...
    num2str(ePSpec,'%8.3e')],'fontsize', 12)
  text(1.25,0.125,['max error fd = ' ...
    num2str(eFDiff,'%8.3e')],'fontsize', 12)
  legend('pseudo-spectral','finite difference','exact')
  
  subplot(2,1,2);
  eigenPSpec(1) = 0.0;
  eigenFDiff(1) = 0.0;
  plot(0:N-1,h*h*eigenPSpec,'-s',0:N-1,h*h*eigenFDiff,'-o')
  grid on, xlabel 'wave number, k', ...
    ylabel 'normalized eigenvalues'
  title(['Normalized Eigenvalues of the Pseudo-Spectral' ...
    ' and Finite Difference Operators']);
  axis([0,N,0,pi*pi])
  set(gca,'xTick',0:N/8:N)
  legend('pseudo-spectral','finite difference')
  
  s1 = ['000' num2str(N)];
  s2 = s1((length(s1)-3):length(s1));
  s3 = ['OUT/periodicPoisson', s2, '.pdf'];
  exportgraphics(hf, s3)

end