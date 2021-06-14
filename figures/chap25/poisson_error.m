e_ps = zeros(10,1);
e_fd = zeros(10,1);

m = [8,12,16,24,32,48,64,96,128,192,256,384];

for k = 1:12
    [e_ps(k),e_fd(k)] = poisson_1d(2.0,m(k),0);
end

comp = m.^(-2);

clf
hff=figure
clf
loglog(m,e_ps,'s',m,e_fd,'o',m,comp,'k')
grid on, xlabel "grid size, $M$", ylabel "error"
title(["Log--Log Plot of Pseudo-Spectral and Finite" ...
    " Difference Errors"]);
axis([5 600 1e-16 1])
legend('pseudo--spectral','finite difference',"$M^{-2}$")
print(hff,"OUT/poisson_error",'-dpdflatex')