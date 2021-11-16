clear all
syms x k j

k = -10:1:10;

ak = int(x*exp(-2*pi*sqrt(-1)*k*x),x,0,0.5)...
    + int((1-x)*exp(-2*pi*sqrt(-1)*k*x),x,0.5,1);

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

figure(1)
clf

plot(z,real(f_proj_3),z,real(f_proj_5),z,real(f_proj_9))

hold on

y = (sawtooth(z*2*pi,1/2)+1)/4;

plot(z,y)

xlabel('$x$','Interpreter','latex');
ylabel('$\mathcal{S}_n[f](x)$','Interpreter','latex');
title('Fourier Projections of a Triangle Wave');
legend('$n = 3$', '$n = 5$', '$n = 9$','Interpreter','latex')
printstr = 'TriangleWave.pdf';
exportgraphics(gca, printstr)

figure(2)
clf

plot(z,y-real(f_proj_3),z,y-real(f_proj_5),z,y-real(f_proj_9))

xlabel('$x$','Interpreter','latex');
ylabel('$f(x)-\mathcal{S}_n[f](x)$','Interpreter','latex');
title('Errors in the Fourier Projections of a Triangle Wave');
legend('$n = 3$', '$n = 5$', '$n = 9$','Interpreter','latex')
printstr = 'TriangleWaveError.pdf';
exportgraphics(gca, printstr)


