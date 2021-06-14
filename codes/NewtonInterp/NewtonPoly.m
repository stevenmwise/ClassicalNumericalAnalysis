x = [0.1,0.3,0.4,0.7,0.9];
y = [5,-1,0,2,2];
a = [5,-30,400/3,-2125/9,5125/18];
z = 0:0.0001:1.0;

i3 = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1)).*(z-x(2)) ...
    + a(4)*(z-x(1)).*(z-x(2)).*(z-x(3)) ;

i4 = i3 + a(5)*(z-x(1)).*(z-x(2)).*(z-x(3)).*(z-x(4));

plot(x,y,'o',z,i3,z,i4)

xlabel('$x$','Interpreter','latex');
title('Third and Fourth-Order Interpolating Polynomials');
legend('$f(x_i)$', '$\mathcal{I}_{X_3}[f](x)$', '$\mathcal{I}_{X_4}[f](x)$','Interpreter','latex')
printstr = strcat('NewtonInterp.pdf');
exportgraphics(gca, printstr)
    
