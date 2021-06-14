clear all;
close all;
f = @(x) 1./(1 + 25*x.^2);

figure(1);
clf;

x = linspace(-1,1,500);
y_true = f(x);
plot(x,y_true,'r','linewidth',2);
hold on;

N = 4;
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_04 = plot(x,y_fit,'g','linewidth',2);


N = 8;
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_08 = plot(x,y_fit,'c','linewidth',2);

N = 16;
xdata = linspace(-1,1,N+1)';
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_16 = plot(x,y_fit,'k','linewidth',2);

plot(xdata,ydata,'k.','markersize',30);

lgd = legend('$f(x)$','$p_4(x)$','$p_8(x)$','$p_{16}(x)$', ...
    'interpolation points','Location','south')
lgd.FontSize = 14;
lgd.Interpreter = 'latex';

axis([-1 1 -1 1])

figure(2);
clf;

x = linspace(-1,1,500);
y_true = f(x);
plot(x,y_true,'r','linewidth',2);
hold on;

N = 4;
idx = 1:N+1;
xdata = cos((idx-0.5)*pi/(N+1));
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_04 = plot(x,y_fit,'g','linewidth',2);
plot(xdata,ydata,'g.','markersize',20);


N = 8;
idx = 1:N+1;
xdata = cos((idx-0.5)*pi/(N+1));
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_08 = plot(x,y_fit,'c','linewidth',2);
plot(xdata,ydata,'c.','markersize',20);

N = 16;
idx = 1:N+1;
xdata = cos((idx-0.5)*pi/(N+1));
ydata = f(xdata);
p = polyfit(xdata,ydata,N);
y_fit = polyval(p,x);
poly_16 = plot(x,y_fit,'k','linewidth',2);
plot(xdata,ydata,'k.','markersize',20);

lgd = legend('$f(x)$','$p_4(x)$','interpolation points','$p_8(x)$', ...
    'interpolation points','$p_{16}(x)$','interpolation points','Location','south')
lgd.FontSize = 14;
lgd.Interpreter = 'latex';

axis([-1 1 -1 1])

figure(3);
clf;
x = linspace(-1,1,500);
N = 8;
idx = 1:N+1;
xdata_cheb = cos((idx-0.5)*pi/(N+1));
xdata_unif = linspace(-1,1,N+1)';

for k = 1:500
    y_cheb(k) = (x(k)-xdata_cheb(1))*(x(k)-xdata_cheb(2))* ...
                (x(k)-xdata_cheb(3))*(x(k)-xdata_cheb(4))* ...
                (x(k)-xdata_cheb(5))*(x(k)-xdata_cheb(6))* ...
                (x(k)-xdata_cheb(7))*(x(k)-xdata_cheb(8))*(x(k)-xdata_cheb(9));
    
    y_unif(k) = (x(k)-xdata_unif(1))*(x(k)-xdata_unif(2))* ...
                (x(k)-xdata_unif(3))*(x(k)-xdata_unif(4))* ...
                (x(k)-xdata_unif(5))*(x(k)-xdata_unif(6))* ...
                (x(k)-xdata_unif(7))*(x(k)-xdata_unif(8))*(x(k)-xdata_unif(9));
    y_bndd(k) = 2^(-N);
end

plot(x,y_cheb,'linewidth',1.0)
hold on
plot(x,y_unif,'linewidth',1.0)
plot(x, y_bndd,'k-.','linewidth',2)
plot(x,-y_bndd,'k:','linewidth',2)

lgd = legend('$\omega_9(x)$ (Chebyshev)','$\omega_9(x)$ (Uniform)',...
             '$2^{-8}$','$-2^{-8}$')
lgd.FontSize = 14;
lgd.Interpreter = 'latex';


