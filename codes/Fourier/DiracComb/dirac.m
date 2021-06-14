syms x
c = int((1+cos(2*pi*x))^{20}, x, 0, 1)
z = -2:0.001:2;
y = (1+cos(2*pi.*z)).^20;
plot(z,y/c,'b-','LineWidth',1.5)
%
xlabel('$x$','Interpreter','latex');
ylabel('$P_{20}(x)/\int_0^1 P_{20}(x) \, d x$','Interpreter','latex');
title('An Approximation of the Periodic Delta Function');
printstr = strcat('DiracComb20.pdf');
exportgraphics(gca, printstr)