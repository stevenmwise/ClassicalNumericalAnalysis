syms x
%
n = 41
%
format long
c = zeros(n,1);
%
for k = 1:2:n
    c(k) = (2*k+1)/2*int(legendreP(k,x),x,0,1);
end
%
z = -1:0.001:1;
lz = length(z);
%
f = ones(1,lz);
f = 0.5*f;
g = zeros(1,lz);
g(1001) = 0.5;
g(1002:lz) = 1.0;
%
for k = 1:2:n
    f = f+c(k)*legendreP(k,z);
end
%
plot(z(1:1000),g(1:1000),'r-','LineWidth',1.5)
hold on
plot(z(1002:lz),g(1002:lz),'r-','LineWidth',1.5)
hold on
plot(z,f,'b-','LineWidth',1.5)
%
xlabel('x');
title(['Fourier-Legendre Expansion of Degree ', num2str(n)]);
axis([-1.0,1.0,-0.1,1.1])
printstr = strcat('FourierLegendre',num2str(n),'.pdf');
exportgraphics(gca, printstr)