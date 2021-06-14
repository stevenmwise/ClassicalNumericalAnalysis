

quad_options("absolute tolerance",eps);
c = quad(@(x) power( 1 + cos(2.*pi.*x), 20), 0, 1);

z = -2:0.001:2;

y = (1+cos(2*pi.*z)).^20;
hf = figure(1);
clf;
plot(z,y./c,'b-','LineWidth',1.5)
%
xlabel('$x$');
ylabel('$P_{20}(x)/\int_0^1 P_{20}(x) \diff x$');
title('An Approximation of the Periodic Delta Function');
grid on;
printstr = strcat('OUT/DiracComb20');
%  exportgraphics(gca, printstr)
print(hf,printstr,'-dpdflatex')
