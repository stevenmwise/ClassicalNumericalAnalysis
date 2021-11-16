% A procedure to plot the boundary of the linear stability 
% domains for the BDF2, BDF3, and BDF4 methods:
t = 0:0.01:2*pi;
%
% BDF2
x2 = 1.5-2.0*cos(t)+0.5*cos(2.0*t);
y2 = 2.0*sin(t)-0.5*sin(2.0*t);
%
% BDF3
x3 = 11/6-3*cos(t)+3/2*cos(2*t)-1/3*cos(3*t);
y3 = 3*sin(t)-3/2*sin(2*t)+1/3*sin(3*t);
%
% BDF4
x4 = 25/12-4*cos(t)+3*cos(2*t)-4/3*cos(3*t)+1/4*cos(4*t);
y4 = 4*sin(t)-3*sin(2*t)+4/3*sin(3*t)-1/4*sin(4*t);

plot(x2,y2,x3,y3,x4,y4)
grid on
