clear all
syms x k j

k = -20:1:20;

%aj = int(1.0*exp(-2*pi*sqrt(-1)*j*x),x,0,0.5)...
    + int(-1.0*exp(-2*pi*sqrt(-1)*j*x),x,0.5,1);

%simplify(aj)


ak = int(1.0*exp(-2*pi*sqrt(-1)*k*x),x,0,0.5)...
    + int(-1.0*exp(-2*pi*sqrt(-1)*k*x),x,0.5,1);

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

figure(1)
clf

plot(z,real(f_proj_15),z,real(f_proj_17),z,real(f_proj_19))

hold on

y = square(2.0*pi*z);

plot(z,y)

xlabel('$x$','Interpreter','latex');
ylabel('$\mathcal{S}_n[f](x)$','Interpreter','latex');
title('Fourier Approximations of a Square Wave','Interpreter','latex');
legend('$n = 15$', '$n = 17$', '$n = 19$','Interpreter','latex')
printstr = strcat('SquareWave.pdf');
exportgraphics(gca, printstr)

figure(2)
clf

plot(z,y-real(f_proj_15),z,y-real(f_proj_17),z,y-real(f_proj_19))

xlabel('$x$','Interpreter','latex');
ylabel('$f(x)-\mathcal{S}_n[f](x)$','Interpreter','latex');
title('Errors in the Fourier Approximations of a Square Wave','Interpreter','latex');
legend('$n = 15$', '$n = 17$', '$n = 19$','Interpreter','latex')
printstr = strcat('SquareWaveError.pdf');
exportgraphics(gca, printstr)


