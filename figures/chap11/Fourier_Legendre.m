pkg load miscellaneous
n = 41
%
c = zeros(n,1);
%
quad_options("absolute tolerance",eps);
for k = 1:2:n
  c(k) = 0.5*(2*k+1)*quad(@(x) legendrepoly( k, x ), 0, 1);
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
    f = f+c(k).*legendrepoly(k,z);
end
%
hf = figure(1);
clf;

plot(z(1:1000),g(1:1000),'r-','LineWidth',1.5)
hold on
plot(z(1002:lz),g(1002:lz),'r-','LineWidth',1.5)
hold on
plot(z,f,'b-','LineWidth',1.5)
%
xlabel('$x$');
title(['Fourier--Legendre Expansion of Degree ', num2str(n)]);
axis([-1.0,1.0,-0.1,1.1])
grid on;
printstr = strcat('OUT/FourierLegendre',num2str(n));
print(hf,printstr,'-dpdflatex')
