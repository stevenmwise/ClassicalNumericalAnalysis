clear all
pkg load signal

ak = zeros(21,1);
quad_options("absolute tolerance",eps);
for k = -10:1:10;
  ak(k+11) = quad( @(x) 0.25*( sawtooth(x*2*pi,0.5) + 1. )*cos(-2*pi*k*x), 0, 1 ) + sqrt(-1)*( quad( @(x) 0.25*( sawtooth(x*2*pi,0.5) + 1. )*sin(-2*pi*k*x), 0, 1 ) );
end
    
z = 0:0.001:2;

lz = length(z);

f_proj_9 = zeros(1,lz);

for kk = -9:1:9
    f_proj_9 = f_proj_9+ak(kk+11)*exp(2*pi*sqrt(-1)*kk*z);
end

f_proj_5 = zeros(1,lz);

for kk = -5:1:5
    f_proj_5 = f_proj_5+ak(kk+11)*exp(2*pi*sqrt(-1)*kk*z);
end

f_proj_3 = zeros(1,lz);

for kk = -3:1:3
    f_proj_3 = f_proj_3+ak(kk+11)*exp(2*pi*sqrt(-1)*kk*z);
end

hf = figure(1)
clf

%  plot(z,real(f_proj_3),z,real(f_proj_5),z,real(f_proj_9))
plot(z,f_proj_3,z,f_proj_5,z,f_proj_9)

hold on

y = (sawtooth(z*2*pi,1/2)+1)/4;

plot(z,y)

xlabel('$x$');
ylabel('$\mathcal{S}_n[f](x)$');
title('Fourier Projections of a Triangle Wave');
legend('$n = 3$', '$n = 5$', '$n = 9$')
grid on;
printstr = strcat('OUT/TriangleWave');
%  exportgraphics(gca, printstr)
print(hf,printstr,'-dpdflatex')


hff=figure(2)
clf

%  plot(z,y-real(f_proj_3),z,y-real(f_proj_5),z,y-real(f_proj_9))
plot(z,y-f_proj_3,z,y-f_proj_5,z,y-f_proj_9)
xlabel('$x$');
ylabel('$f(x)-\mathcal{S}_n[f](x)$');
title('Errors in the Fourier Projections of a Triangle Wave');
legend('$n = 3$', '$n = 5$', '$n = 9$')
grid on;
printstr = strcat('OUT/TriangleWaveError');
%  exportgraphics(gca, printstr)
print(hff,printstr,'-dpdflatex')


