clf
clear

xf = -0.1:0.005:2.1;
yf = exp(sin(2*pi*xf)-0.1*cos(2*pi*xf)+0.2*cos(4*pi*xf)+0.2*sin(6*pi*xf));

pn =0;

for n = [6,8,12,16]
    pn = pn+1;

    h = 1/n;

    xc = 0.0:h:1-h;
    yc =  exp(sin(2*pi*xc)-0.1*cos(2*pi*xc)+0.2*cos(4*pi*xc)+0.2*sin(6*pi*xc));

    n = length(yc);
    m = floor((n+1)/2);
    z = fft(yc)/n;

    a0 = z(1); 
    ak(1:m-1) = 2*real(z(2:m));
    ak(m) = z(m+1);
    bk(1:m-1) = -2*imag(z(2:m));

    k = 1:length(ak)-1;

    py = a0 + ak(1:m-1)*cos(2*pi*k'*xf)...
        + bk(1:m-1)*sin(2*pi*k'*xf) ...
        + ak(m)*cos(2*pi*m*xf); 
    
    hf = figure(pn)
    clf
    
    plot(xf,yf,'b-',xc,yc,'ro',xf,py,'r--','LineWidth',1.25)
    axis([-0.1,2.1,0.2,2.6])
    grid on
    xlabel('$x$');
    title(['Trigonometric Interpolation, $n=$',num2str(n)]);
    legend('$f(x)$', '$f(x_j)$', 'interpolant')
    printstr = strcat('OUT/TrigInterp',num2str(n),'.pdf')
    
    print( hf, printstr, '-dpdflatex' );
    
    printstrGray = strcat('OUT/TrigInterp',num2str(n),'Gray.pdf');
    print( hf, printstrGray, '-dpdflatex', '-mono')
end
