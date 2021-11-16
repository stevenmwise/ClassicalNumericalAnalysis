clear all;
clc;
close all;


degree = 20;
step = 1e-4;
pert = 1e-9;

hf = figure();
plot(zeros(degree,1),'gr *;roots;',"linewidth",1);
xlabel("$t$");
ylabel("$p(t)$");
grid on
hold on;
title( "The Wilkinson polynomial and its perturbation" );

x=(1-step):step:(degree+step);

yexact = WilkinsonFact(x);
legend("location","southeast");
plot(x,yexact,'linewidth',1,";$p_W$;");
hold on;

ynoise = yexact - pert.*(x.^19);
plot(x,ynoise,"linewidth",1,";$\\tilde{p}_W$;");
hold on

print( hf, "OUT/WilkinsonPlot", "-dpdflatex");

hfgray = hf;
print( hfgray, "OUT/WilkinsonPlotGray", "-dpdflatex", "-mono");
