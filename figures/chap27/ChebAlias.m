N = 8;
M = 2001;
x = [-1,-cos(pi*(1:N)./N)];
xp = linspace(-1,1,M);

T21 = cos(21*acos(xp));
T11 = cos(11*acos(xp));
T5  = cos( 5*acos(xp));
T5g = cos( 5*acos(x ));

XCGL = zeros(1,length(x));

hf = figure(1);
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

xlabel('$x$');
ylabel('$T_5(x),\ T_{11}(x),\ T_{21}(x)$');
title('Aliasing of Chebyshev Polynomials $T_5,\ T_{11},\ T_{21}$, $N=8$');
axis([-1.0,1.0,-1.1,1.1])
grid on
printstr = strcat('OUT/ChebyshevAlias',num2str(N));
printstrGray = strcat('OUT/ChebyshevAlias',num2str(N),'Gray');
%exportgraphics(gca, printstr)
print(hf, printstr, '-dpdflatex');
print(hf, printstrGray, '-dpdflatex','-mono');
