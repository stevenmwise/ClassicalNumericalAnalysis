  x0 = 0.0:0.01:1.0;
y1 = 4.0.*sin(pi*x0).*sin(pi*x0);
y2(1:101) = 4.0;
y3(1:101) = 0.0;
h = 1/10;
x1 = h:h:1;
y4 = 4.0.*sin(pi*x1).*sin(pi*x1);
plot(x0,y1,'k-',x0,y2,'k--',x0,y3,'k--',x1,y4,'o')
axis([0,1,-1,5])
print('eigenvalues.eps','-deps')