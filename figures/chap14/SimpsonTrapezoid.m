function [] = quad() 

close all
clear all
clf

m = 4;
zz = linspace(0,pi,1001);
yy = f(zz);
exact = -12.070346316389633;

hf = figure(1);
clf;
z = linspace(0,pi,m+1);
format long
trap = (0.5*f(z(1))+0.5*f(z(m+1))+sum(f(z(2:m))))*pi/m
trapcorr = trap-(1/12)*(f1(pi)-f1(0))*(pi/m)^2
format longe
etrap = exact-trap
etrapcorr = exact-trapcorr

y = f(z);


h = plot(zz,yy,'k',z,y,'--r',z,y,'ko')
grid on;
set(h(1),'linewidth',1);
set(h(2),'linewidth',2);
set(h(3),'MarkerFaceColor', 'k');

axis([0,pi,-24,4])
xlabel('$x$');
title(['Composite Trapezoidal Rule: $m = $ ', num2str(m)]);
printstr = strcat('OUT/Trapezoidal',num2str(m));
printstrGray = strcat('OUT/Trapezoidal',num2str(m),'Gray');
% exportgraphics(gca, printstr)
print(hf,printstr,'-dpdflatex')
print(hf,printstrGray,'-dpdflatex','-mono')

z = linspace(0,pi,2*m+1);

format long
simp = (f(z(1))+f(z(2*m+1))+2.0*sum(f(z(3:2:2*m-1)))+4.0*sum(f(z(2:2:2*m))))*pi/(6*m)
simpcorr = simp-(1/180)*(f3(pi)-f3(0))*(pi/(2*m))^4
format longe
esimp = exact-simp
esimpcorr = exact-simpcorr

hff= figure(2)
clf
for k = 1:m
    zdata = z(1+(k-1)*2:3+(k-1)*2);
    ydata = f(zdata);
    p = polyfit(zdata,ydata,2);
    zzz = linspace(z(1+(k-1)*2),z(3+(k-1)*2),500);
    y_fit = polyval(p,zzz);
    h = plot(zzz,y_fit,'--r',zdata,ydata,'ko')
    set(h(1),'linewidth',2);
    set(h(2),'linewidth',1);
    hold on
end
plot(zz,yy,'k')
hold on
z = linspace(0,pi,m+1);
y = f(z);
plot(z,y,'ko','Markerface','k')
grid on;

axis([0,pi,-24,4])
xlabel('$x$');
title(['Composite Simpson''s Rule: $m=$ ', num2str(m)]);
printstr = strcat('OUT/Simpson',num2str(m));
printstrGray = strcat('OUT/Simpson',num2str(m),'Gray');
% exportgraphics(gca, printstr)
print(hff,printstr,'-dpdflatex')
print(hff,printstrGray,'-dpdflatex','-mono')

hold off
end


function [z] = f(x)
z = exp(x).*cos(x);
end

function [z] = f1(x)
z = exp(x).*(cos(x)-sin(x));
end

function [z] = f3(x)
z = -2.0*exp(x).*(cos(x)+sin(x));
end
