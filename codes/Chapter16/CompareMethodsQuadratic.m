% This code defines a strongly convex and globally Lipschitz 
% smooth quadratic objective and uses it to compare the 
% performance of PSD and AGD for its minimization:
clear all
clc

% We begin by defining the space dimension and the parameters 
% that are used to define the objective function:
dim = 100000; % the dimension: the number of variables
mu = 1.; % the strong convexity constant
L = 1000.; % the global Lipschitz constant

% Since the objective is quadratic, part of it is defined as
%
% v'*AA*v
%
% where AA is a tridiagonal matrix that we now define.
mult = 0.25*( L - mu );
ii = zeros(dim + 2*( dim - 1 ), 1 );
jj = ii;
vv = ii;
vvEuler = ii;
for i = 1:dim
  ii(i) = i;
  jj(i) = i;
  vv(i) = 2;
  vvEuler(i) = mult*vv(i) + mu;
  if i<dim
    ii(i+dim) = i;
    jj(i+dim) = i+1;
    vv(i+dim) = -1;
    vvEuler(i+dim) = mult*vv(i+dim);
  end
  
  if i>1
    ii(i+2*(dim-1)) = i;
    jj(i+2*(dim-1)) = i-1;
    vv(i+2*(dim-1)) = -1;
    vvEuler(i+2*(dim-1)) = mult*vv(i+2*(dim-1));
  end
end

AA = sparse(ii,jj, vv );
AAEuler = sparse( ii, jj, vvEuler );

e1 = zeros(dim,1);
e1(1,1) = 1;

% We define two anonymous function handles to describe the 
% objective and its derivative:
E = @(u) Objective( u, mu, L, AA, e1 );
Ep = @(u) ObjectivePrime( u, mu, L, AA, e1 );

% The minimum of the objective can be found by solving the 
% Euler equations. Since the objective is quadratic these turn 
% out to be a linear system of equations:
uexact = AA*e1;
uexact = AAEuler\uexact;
uexact = 0.25*( L - mu )*uexact;
Eexact = E( uexact );

% We now can start with the optimization schemes. Maximum 
% number of iterations and initial guess:
MaxIts = 1000; % Maximum number of iterations
uinit = zeros( dim, 1 ); % Initial guess

% We run the methods
[energyPSD, energyAGD, upsd, uagd] = RunMethods( E, Ep, mu, ...
  L, MaxIts, uinit );
energyPSD = energyPSD - Eexact;
energyAGD = energyAGD - Eexact;

% We plot to compare
its = 0:MaxIts;
hf = figure();
clf;
plot( its, log10( energyPSD ), its, log10( energyAGD ) );
grid on;
xlabel('k');
ylabel('log_{10}( E(u_k) - E(u) )');
title('Comparison of PSD and AGD for a quadratic objective');
legend('PSD', 'AGD');
exportgraphics( gca, 'ComparisonQuadratic.pdf');

% We now compare the rate of convergence of AGD with the 
% theoretically optimal one:
hff = figure(2);
clf;
kappa = sqrt( L/mu );
rate =  2.*log10( 1 - 1/kappa )*its;

plot( its, log10( energyAGD ), its, rate );
grid on;
xlabel('k');
ylabel('log_{10}( E(u_k) - E(u) )');
legend('AGD','optimal rate')
title('Error of AGD vs. Optimal rate');
exportgraphics( gca,'ComparisonAGDandOptimal.pdf');

%%%%%%%%%%%%%%%%%%%%%%%

function E = Objective( u, mu, L, AA, e1 )
% This function defines a quadratic objective
%
% Input
%   u : the argument of the Objective function
%   mu : the strong convexity constant of the objective
%   L : the Lipschitz smoothness constnat of the objective
%   AA : an SPD matrix used to define the objective
%   e1 : A fixed vector to define the objective
%
% Output
%   E : the value of the objective function
%
  E = 0.125*( L - mu )*( ( u - e1 )'*AA*( u - e1 ) ) ...
    + 0.5*mu*( norm(u)^2 );
end

function Ep = ObjectivePrime( u, mu, L, AA, e1 )
% This function defines the derivative of a quadratic objective
%
% Input
%   u : the argument of the Objective function
%   mu : the strong convexity constant of the objective
%   L : the Lipschitz smoothness constnat of the objective
%   AA : an SPD matrix used to define the objective
%   e1 : A fixed vector to define the objective
%
% Output
%   Ep : the value of the derivative of the objective function
%
  Ep = 0.25*(L-mu)*AA*(u - e1) + mu*u;
end
