%% A very simple procedure to plot the boundary of the linear stability domain of BDF2
t = 0:0.01:2*pi;

x = 1.5-2.0*cos(t)+0.5*cos(2.0*t);
y = 2.0*sin(t)-0.5*sin(2.0*t);

plot(x,y)
