clear all;
clc;
close all;

hf = figure();
%axis equal
%axis([-1,7,-4,4])
xlabel("$\\Re(z)$");
ylabel("$\\Im(z)$");
grid on
hold on;
%title( "The boundary of the linear stability domain $\\mathcal{D}_{\\textrm{BDF2}}$ for the BDF2 method." );

t = 0:0.01:1.0;
t = 2.*pi.*t;
x = 1.5 - 2.*cos(t)+0.5.*cos(2.*t);
y = 2*sin(t)-0.5*sin(2.*t);
plot(x,y)
hold on;

print( hf, "BDF2BdryLocus.pdf", "-dpdflatex");

