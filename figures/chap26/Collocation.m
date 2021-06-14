figure(1)
clf
figure(2)
clf
figure(1)

N = [8,12,16,20,24,32,40,48,64,80];

for i = 1:length(N)
    [ecoll(i),efd(i)] = BVPCollVsFD(N(i));
end

comp = 4*N.^(-2);

htr= figure(2);
loglog(N,ecoll,'ro',N,efd,'bs',N,comp,'k-')

grid on
xlabel('$N$')
ylabel('max norm error at the nodes')
title(['Log-Log Plots of Collocation and Finite Difference Errors']);
axis([7 90 2e-15 2])
legend('Chebyshev Collocation','Finite Difference', '$4N^{-2}$','Location','east')

% exportgraphics(gca, 'CollocationErrors.pdf')
print(htr,'OUT/CollocationErrors','-dpdflatex');
