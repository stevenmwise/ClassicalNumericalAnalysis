clear all;
clc;

z =linspace(0.0,1.0,8001);

y = ones(size(z));
w = zeros(size(z));
y = (z-0*0.1250).*(z-1*0.1250).*(z-2*0.1250).*(z-3*0.1250).*(z-4*0.1250).* ...
    (z-5*0.1250).*(z-6*0.1250).*(z-7*0.1250).*(z-8*0.1250);
%
% Output plots:
hf = figure(1);
clf;

H=area(z(0001:1000),y(0001:1000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(1001:2000),y(1001:2000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(2001:3000),y(2001:3000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(3001:4000),y(3001:4000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(4001:5000),y(4001:5000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(5001:6000),y(5001:6000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(6001:7000),y(6001:7000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(7001:8001),y(7001:8001));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

plot(z,y,'k',z,w,'k')

ylabel('$\omega_9(x)$')
xlabel('$x$');
title('The Nodal Polynomial $\omega_9$ with Equally Spaced Nodes');

printstr = strcat('OUT/NodalPolynomial');
%  exportgraphics(gca, printstr)
print(hf,printstr,'-dpdflatex')

hff = figure(2);
clf;
inty = cumsum(y)*(1/8000);

H=area(z(0001:1000),inty(0001:1000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(1001:2000),inty(1001:2000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(2001:3000),inty(2001:3000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(3001:4000),inty(3001:4000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(4001:5000),inty(4001:5000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(5001:6000),inty(5001:6000));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on

H=area(z(6001:7000),inty(6001:7000));
set(H(1),'FaceColor',[0.85 0.85 0.85]);
hold on

H=area(z(7001:8001),inty(7001:8001));
set(H(1),'FaceColor',[0.5 0.5 0.5]);
hold on


plot(z,inty,'k')

ylabel('$\Omega_9(x)$')
xlabel('$x$');
title('$\Omega_9(x) = \int_0^x\omega_9(t)\, \mbox{d} t$');

printstr = strcat('OUT/IntNodalPolynomial');
%  exportgraphics(gca, printstr)
print(hff,printstr,'-dpdflatex')

        
