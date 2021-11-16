N = 8;
M = 2001;
x = [-1,-cos(pi*(1:N)./N)];
xp = linspace(-1,1,M);

T21 = cos(21*acos(xp));
T11 = cos(11*acos(xp));
T5  = cos( 5*acos(xp));
T5g = cos( 5*acos(x ));

XCGL = zeros(1,length(x));

clf
plot(xp,T5 ,'r','LineWidth',1.5)
hold on
plot(xp,T11,'b','LineWidth',1.5)
hold on
plot(xp,T21,'g','LineWidth',1.5)
hold on
plot(x ,T5g,'ko','LineWidth',1.5)
hold on
plot(x,XCGL,'k-o','LineWidth',1.5)
hold off

xlabel('$x$','Interpreter','latex');
ylabel('$T_5(x),\ T_{11}(x),\ T_{21}(x)$','Interpreter','latex');
title('Aliasing of Chebyshev Polynomials $T_5,\ T_{11},\ T_{21}$, $N=8$', ...
    'Interpreter','latex');
axis([-1.0,1.0,-1.1,1.1])
printstr = strcat('ChebyshevAlias',num2str(N),'.pdf');
exportgraphics(gca, printstr)