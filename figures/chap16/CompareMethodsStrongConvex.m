% This code defines an objective that not quadratic, it is strongly convex,
% and locally but not globally Lipschitz smooth. This objective is used to
% compare the performance of PSD and AGD
clear all
clc

% We begin by defining the space dimension and the parameters that are used
% to define the objective function
dim = 1000000; % the dimension: the number of variables
p = 4; % this is used to define the norm
mu = 1; % the strong convexity constant
% we do not know exactly the (local) Lipschitz smoothness constant. This is
% an estimate
L = 400;

%%%%%%%%%%%%%%%%%%%%%%%

function E = Objective( u, p )
% This defines a non quadratic, strongly convex, and locally but not
% globally Lipshitz smooth objective.
%
% Input
%   u : the argument of the Objective function
%   p : the order of the norm used to make this non quadratic
%
% Output
%   E : the value of the objective function
%
  E = ( norm( u, p )^p )/p + 0.5*norm( u )^2;
end
  
function Ep = ObjectivePrime( u, p, dim )
% This defines the derivative of the nonquadratic objective
%
% Input
%   u(1,dim) : the argument of the Objective function
%   p : the order of the norm used to make this non quadratic
%   dim : the dimension of u
%
% Output
%   Ep : the value of the derivative of the objective function
%
  Ep = u + ( abs( u ).^(p-2) ).*u;
end

% We define two anonymous function handles to describe the objective and
% its derivative
E = @(u) Objective( u, p );
Ep = @(u) ObjectivePrime( u, p, dim );

% We know the exact minimizer
uexact = zeros( dim, 1 );
Eexact = E( uexact );

% We now can start with the optimization schemes. Maximum number of iterations and initial guess
MaxIts = 1000; % Maximum number of iterations
uinit = 10*ones(dim,1); % Initial guess

% We run the methods
[errorPSD, errorAGD, upsd, uagd] = RunMethods( E, Ep, mu, L, MaxIts, uinit );
errorPSD = errorPSD - Eexact;
errorAGD = errorAGD - Eexact;

% We plot to compare
its = 0:MaxIts;
hf = figure();
clf;
plot( its, log10( errorPSD ), its, log10( errorAGD ) );
grid on;
xlabel( "$k$" );
ylabel( "$\\log_{10}( E(u_k) - E(u) )$" );
title( "Comparison of PSD and AGD for a strongly convex objective" );
legend("PSD", "AGD" );

print( hf, ["OUT/ComparisonStrongConvex", num2str(dim), ".pdf"], "-dpdflatex");

print( hf, ["OUT/ComparisonStrongConvex", num2str(dim), "Gray.pdf"], "-dpdflatex", "-mono");
