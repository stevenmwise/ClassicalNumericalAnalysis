% Hermite interpolation 
 
n = 6;
x = zeros(1,n+1);
y = zeros(1,n+1);
q = zeros(2*n+2,2*n+2);
 
x(1) =  0.00; q(01,1) = -1.0; q(02,1) = 0.0;
x(2) =  0.20; q(03,1) =  1.0; q(04,2) = 0.0;
x(3) =  0.50; q(05,1) =  1.0; q(06,2) = 0.0;
x(4) =  1.00; q(07,1) = -1.0; q(08,2) = 0.0;
x(5) =  1.50; q(09,1) = -1.0; q(10,2) = 0.0;
x(6) =  2.00; q(11,1) =  1.0; q(12,2) = 0.0;
x(7) =  2.20; q(13,1) = -1.0; q(14,2) = 0.0;
%
y(1:n+1) = q(1:2:2*n+1,1);
xo = [x(1),x(3),x(4),x(6),x(7)];
yo = [y(1),y(3),y(4),y(6),y(7)];
  
 z = zeros(1,2*n+2);
 for i = 0:n
    z(2*i+1) = x(i+1);
    z(2*i+2) = x(i+1);
    q(2*i+2,1) = q(2*i+1,1);
    if i ~= 0 
       q(2*i+1,2) = (q(2*i+1,1)-q(2*i,1))/(z(2*i+1)-z(2*i));
    end
 end

 k = 2*n+1;
 for i = 2:k
    for j = 2:i
       q(i+1,j+1) = (q(i+1,j)-q(i,j))/(z(i+1)-z(i-j+1));
    end
 end

a = min(x)-0.05; b = max(x)+0.05;
M = 500;
h = (b-a)/M;
xx = a:h:b;
yy = zeros(1,length(xx));
for ell = 1:M+1
    s = q(k+1,k+1)*(xx(ell)-z(k));
    for i = 2:k
        j = k-i+1;
        s = (s+q(j+1,j+1))*(xx(ell)-z(j));
    end
    s = s + q(1,1);
    yy(ell) = s;
end
clf;
zo = zeros(1,length(xo));
plot(xx,yy,'k-',xo,yo,'ko','MarkerFaceColor','k')
hold on;
plot(xo,zo,'ko')
axis([a b -1 1]);
grid on;
xi(1) = 0.1; xi(2) = 0.75; xi(3) = 1.75; xi(4) = 2.1;
for ell = 1:4
    s = q(k+1,k+1)*(xi(ell)-z(k));
    for i = 2:k
        j = k-i+1;
        s = (s+q(j+1,j+1))*(xi(ell)-z(j));
    end
    s = s + q(1,1);
    yi(ell) = s;
end

zi = zeros(1,length(xi));
plot(xi,yi,'ks','MarkerFaceColor','k')
plot(xi,zi,'ks')
lgd = legend('$f(x)-p(x)$','$f(x_i)-p(x_i)$','$x_i$', ...
  '$f(\xi_i)-p(\xi_i)$','$\xi_i$','Location','north')
lgd.FontSize = 14;
lgd.Interpreter = 'latex';

hold off;
exportgraphics(gca, 'Oscillation.pdf')



