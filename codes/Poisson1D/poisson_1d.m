function [e_ps,e_fd] = poisson_1d(L,m)
%
% This function computes pseudo-spectral and finite difference 
% approximations to the periodic Poisson equation:
%
%   -D_{xx} u = f, where u and f are [0,L]-periodic
%
% using the Fast Fourier Transform algorithm.
%
% Input
%
%   L : the size of the domain
%   m : the number of points in the periodic grid. Generally this should be a positive, even number.
%
% Output
%   e_ps : the max norm of the error in the pseudo-spectral computation.
%   e_fd : the max norm of the error in the finite difference computation.
%
%%%%%%%%%
% h is the grid spacing and x holds the nodal values of the grid.
h = L/m;
x = (1:m)*h;
%
% Fine grid parameters for representing the exact solution.
mf = 256;
hf = L/mf;
xf = (1:mf)*hf;
%
% q is a useful factor that allows us to easily change the 
%   domain size.
q = 2*pi/L;
%
% Wave numbers, k: This definition is really convenient,
%   essentially making k periodic. Note that our array index
%   starts at 1, as usual in MATLAB, but the value of k starts
%   at 0!
k = [0:m/2-1 m/2 -m/2+1:-1];
%
% Right-Hand-Side function, f, forced to be discrete mean-
%   zero:
f = q*q*(cos(q*x).*cos(q*x)-sin(q*x)).*exp(sin(q*x));
f = -(f-h*sum(f)/L);
%
% Exact solution, forced to be discrete mean-zero:
u_exact = exp(sin(q*x));
u_exact_mass = h*sum(u_exact)/L;
u_exact = u_exact-u_exact_mass;
%
% Eigenvalues of the pseudo-spectral derivative operator -D^2:
eigen_ps = q*q*k.*k;
%
% Modified to avoid dividing by zero below:
eigen_ps(1) = 1.0;
%
% Eigenvalues of the finite difference operator -\Delta_h: 
%   Note that we have applied a shift here to be consistent 
%   with the pseudo-spectral case.
eigen_fd = 4.0*sin(q*(x-h)/2).*sin(q*(x-h)/2)/(h*h);
%
% Modified to avoid dividing by zero below:
eigen_fd(1) = 1.0;
%
% Approximations: Note that MATLAB fft's do not respect the
%   fact that our data points are real numbers. So we need to
%   get rid of any imaginary numbers in the results.
u_ps = real((ifft(fft(f)./eigen_ps)));
u_fd = real((ifft(fft(f)./eigen_fd)));
%
% Error computations:
e_ps = max(abs(u_ps-u_exact));
e_fd = max(abs(u_fd-u_exact));
%
% Redefine u_exact for plotting:
u_exact = exp(sin(q*xf));
u_exact = u_exact-u_exact_mass;
%
subplot(2,1,1);
plot(x,u_ps,'s',x,u_fd,'o',xf,u_exact,'k-')
grid on, xlabel x, ylabel 'exact and approximate solutions';
title(['Pseudo-Spectral and Finite Difference ', ...
    'Approximations: m = ', num2str(m)]);
axis([0,L,-1.1,1.7])
set(gca,'xTick',0:L/8:L)
text(1.25,0.425,['max error ps = ' num2str(e_ps,'%8.3e')], ...
    'fontsize', 12)
text(1.25,0.125,['max error fd = ' num2str(e_fd,'%8.3e')], ...
    'fontsize', 12)
legend('pseudo-spectral','finite difference','exact')
%
subplot(2,1,2);
eigen_ps(1) = 0.0;
eigen_fd(1) = 0.0;
plot(0:m-1,h*h*eigen_ps,'-s',0:m-1,h*h*eigen_fd,'-o')
grid on, xlabel 'wave number, k', ...
    ylabel 'normalized eigenvalues'
title(['Normalized Eigenvalues of the Pseudo-Spectral' ...
    ' and Finite Difference Operators']);
axis([0,m,0,pi*pi])
set(gca,'xTick',0:m/8:m)
legend('pseudo-spectral','finite difference')
%
end
