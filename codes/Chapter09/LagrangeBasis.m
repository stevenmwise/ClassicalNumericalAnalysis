function [] = LagrangeBasis(x,namestr)
%
% This function plots the Langrange basis functions in the 
%   interval [0,1] for the interpolation nodes contained
%   in the vector x.
%
% Input:
%
% x is a vector of size n+1 containing the interpolation
%   nodes. These should include 0 and 1 as the first and
%   last nodes and should be strictly increasing.
%
% str is a string that is used in the names of the plot
%   files on output.
%
clf
z= 0.0:0.00005:1.0;
%
n = length(x)-1;
m = length(z);
f = zeros(1,m);
y = eye(n+1);
summ1 = zeros(1,m);
summ2 = zeros(1,m);
%
str = ["-","--",":","-."];
%
maxx = 0;
minn = 0;
figure(1)
for k = 1:n+1
    f(1:m) = polyinterp(x,y(k,1:n+1),z);
    mc = (m-1)/10+1;
    fc(1:mc) = f(1:10:m);
    zc(1:mc) = z(1:10:m);
    plot(x,y(k,1:n+1),'ok',zc,fc(1:mc),'-','LineWidth',1.5)
    hold on 
    summ1(1:m) = summ1(1:m)+f(1:m);
    summ2(1:m) = summ2(1:m)+abs(f(1:m));
    maxx = max(maxx,max(f));
    minn = min(minn,min(f));
end
hold off
%
xlabel('x');
title(['Lagrange Basis Functions of Degree ', num2str(n), ...
    ': ', namestr, ' Nodes']);
axis([0.0,1.0,minn,maxx])
%
% Output plots:
printstr = strcat('Basis0',num2str(n),namestr,'.pdf');
exportgraphics(gca, printstr)
%
figure(2) 
plot(z,summ2(1:m),'-','LineWidth',1.5)
%
xlabel('x');
title(['Lebesgue Function of Degree ', num2str(n), ...
    ': ', namestr, ' Nodes']);
axis([0.0,1.0,1,max(summ2)])
printstr = strcat('Lebesgue0',num2str(n),namestr,'.pdf');
exportgraphics(gca, printstr)
%
end % function

function L = polyinterp(x,y,z)
%
n = length(x);
L = zeros(size(z));
%
for k = 1:n
    w = ones(size(z)); 
    for j = [1:k-1 k+1:n]
        w = (z-x(j))./(x(k)-x(j)).*w; 
    end
    L = L + w*y(k);
end
%
end % function
