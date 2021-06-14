% Collocation

function [errColl,errfd]  = BVPCollVsFD(N)

[D,x] = cheb(N);        
D2 = D^2;              
D2 = D2(2:N,2:N)-eye(N-1); 
f = exp(4*x(2:N))-(exp(4*x(2:N))-sinh(4)*x(2:N)-cosh(4))/16;           
u = D2\f;
u = [0;u;0];
hg = figure(1);
plot(x,u,'r.-','markersize',16)
%
exactColl = (exp(4*x)-sinh(4)*x-cosh(4))/16; 
errColl = norm(u-exactColl,inf);
%
% Finite difference part, for comparison:
%
A = 2.0*eye(N-1);
for i = 1:N-2
    A(i,i+1) = -1.0;
    A(i+1,i) = -1.0;
end
h = 2/N;
A = A+h^2*eye(N-1);
xfd = [-1.0:h:1.0]';
ffd = -h^2*exp(4*xfd(2:N))+h^2*(exp(4*xfd(2:N))-sinh(4)*xfd(2:N)-cosh(4))/16; 
ufd = A\ffd;
ufd = [0;ufd;0];
exactfd = (exp(4*xfd) - sinh(4)*xfd - cosh(4))/16; 
errfd = norm(ufd-exactfd,inf);
hold on
plot(xfd,ufd,'b.-','markersize',16)
xlabel('$x$');
title(['Chebyshev Collocation and Finite Difference Approximations: $N =$ ' num2str(N)]);
legend('Chebyshev Collocation','Finite Difference','Location','southwest')
text(-0.2,-0.5,['max collocation error       = ' num2str(errColl,'%8.3e')], 'fontsize', 12)
text(-0.2,-0.7,['max finite difference error = ' num2str(errfd,  '%8.3e')], 'fontsize', 12)
grid on;
printstr = strcat('OUT/CollVsFD',num2str(N));
% exportgraphics(gca, printstr)
print(hg,printstr,'-dpdflatex');

hold off

end


