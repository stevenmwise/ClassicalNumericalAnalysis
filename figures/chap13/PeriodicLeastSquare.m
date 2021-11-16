% This Matlab script computes the trigonometric interpolant to a
% noisy signal and the associated least squares approximation of 
% lesser degree.
%
clf
clear
%
% The number of noisy samples that are used:
n = 32;
h = 1/n;
%
% n*scale is the number of points that define the uncorrupted and 
% corrupted signals.
scale = 10;
nf = n*scale;
hf = 1/nf;
%
% Clean signal:
xf = -scale*hf:hf:2+scale*hf;
yf = exp(sin(2*pi*xf)-0.1*cos(2*pi*xf)+0.2*cos(4*pi*xf)+0.2*sin(6*pi*xf));
%
% Signal with noise:
seed = 1010;
rand( "seed", seed );
%  rng(seed);
w = 0.075*randn(1,length(xf));
yfw = yf+w;
%
% Get n noisy interpolation/sampe points:
xc = xf( scale+1:scale:nf+scale);
yc = yfw(scale+1:scale:nf+scale);
%
% Get the interpolation coefficients:
m = n/2;
z = fft(yc)/n;
a0 = real(z(1)); 
ak(1:m-1) = 2*real(z(2:m));
ak(m) = real(z(m+1));
bk(1:m-1) = -2*imag(z(2:m));
%
% Construct interpolant:
k = 1:length(ak)-1;
py = a0 + ak(1:m-1)*cos(2*pi*k'*xf)...
        + bk(1:m-1)*sin(2*pi*k'*xf) ...
        + ak(m)*cos(2*pi*m*xf); 
%
% Construct least squares approx:
nls = 8;
mls = nls/2;
als0 = real(z(1)); 
alsk(1:mls-1) = 2*real(z(2:mls));
alsk(mls) = real(z(mls+1));
blsk(1:mls-1) = -2*imag(z(2:mls));
kls = 1:length(alsk)-1;
pyls = als0 + alsk(1:mls-1)*cos(2*pi*kls'*xf)...
          + blsk(1:mls-1)*sin(2*pi*kls'*xf) ...
          + alsk(mls)*cos(2*pi*mls*xf); 

hf = figure();          

plot(xf,yf,'r-','LineWidth',2,xc,yc,'bo','LineWidth',2,xf,py,'b--','LineWidth',2,xf,pyls,'LineWidth',2,'m-.','LineWidth',2)
axis([0.0,1.0,0.2,2.6])
grid on;

xlabel('$x$');
title(['Trigonometric Interpolation, $n=$',num2str(n),' Least Squares $m=$',num2str(nls)]);
legend('$f(x)$', '$f(x_j)+\chi(x_j)$', 'interpolant','least squares')


printstr = strcat("OUT/TrigLeastSquare",num2str(n),".pdf");

print( hf, printstr, '-dpdflatex' );

printstrGray = strcat('OUT/TrigLeastSquare',num2str(n),'Gray.pdf');
print( hf, printstrGray, '-dpdflatex', '-mono')
