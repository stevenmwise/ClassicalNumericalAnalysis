clear
clf

N = 4;
idx = 1:N+1;
xdata = (cos((idx-0.5)*pi/(N+1))+1.0)/2.0;
xdata = flip(xdata);
LagrangeBasis(xdata,'Chebyshev');

N = 8;
idx = 1:N+1;
xdata = (cos((idx-0.5)*pi/(N+1))+1.0)/2.0;
xdata = flip(xdata);
LagrangeBasis(xdata,'Chebyshev');

N = 16;
idx = 1:N+1;
xdata = (cos((idx-0.5)*pi/(N+1))+1.0)/2.0;
xdata = flip(xdata);
LagrangeBasis(xdata,'Chebyshev');