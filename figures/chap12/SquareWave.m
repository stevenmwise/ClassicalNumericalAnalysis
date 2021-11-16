clear all
pkg load signal

ak = zeros(41,1);
quad_options("absolute tolerance",eps);
for k = -20:1:20;
ak(k+21) = ( quad( @(x) cos(-2*pi*k*x), 0, 0.5 ) + quad( @(x) -cos(-2*pi*k*x), 0.5, 1) ) + sqrt(-1)*( quad( @(x) sin(-2*pi*k*x), 0, 0.5 ) + quad( @(x) -sin(-2*pi*k*x), 0.5, 1) );
end
    
z = 0:0.001:2;

lz = length(2.0*pi*z);

f_proj_15 = zeros(1,lz);

for kk = -15:1:15
    f_proj_15 = f_proj_15+ak(kk+21)*exp(2*pi*sqrt(-1)*kk*z);
end

f_proj_17 = zeros(1,lz);

for kk = -17:1:17
    f_proj_17 = f_proj_17+ak(kk+21)*exp(2*pi*sqrt(-1)*kk*z);
end

f_proj_19 = zeros(1,lz);

for kk = -19:1:19
    f_proj_19 = f_proj_19+ak(kk+21)*exp(2*pi*sqrt(-1)*kk*z);
end

hf = figure(1)
clf

plot(z,real(f_proj_15),z,real(f_proj_17),z,real(f_proj_19))

hold on

y = square(2.0*pi*z);

plot(z,y)

xlabel('$x$');
ylabel('$\mathcal{S}_n[f](x)$');
title('Fourier Approximations of a Square Wave');
legend('$n = 15$', '$n = 17$', '$n = 19$')
grid on;
printstr = strcat('OUT/SquareWave');
printstrGray = strcat('OUT/SquareWaveGray');
%  exportgraphics(gca, printstr)
print(hf,printstr,'-dpdflatex')
print(hf,printstrGray,'-dpdflatex','-mono')


hff=figure(2)
clf

plot(z,y-real(f_proj_15),z,y-real(f_proj_17),z,y-real(f_proj_19))

xlabel('$x$');
ylabel('$f(x)-\mathcal{S}_n[f](x)$');
title('Errors in the Fourier Approximations of a Square Wave');
legend('$n = 15$', '$n = 17$', '$n = 19$')
grid on;
printstr = strcat('OUT/SquareWaveError');
printstrGray = strcat('OUT/SquareWaveErrorGray');
%  exportgraphics(gca, printstr)
print(hff,printstr,'-dpdflatex')
print(hff,printstrGray,'-dpdflatex','-mono')


