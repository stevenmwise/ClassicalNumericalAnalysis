function [energyPSD, energyAGD, upsd, uagd] = RunMethods( ...
  E, Ep, mu, L, MaxIts, uinit )
% This function runs both the Steepest Descent with approximate 
% line search (PSD)
%
% u_{k+1} = u_k - alpha*Ep( u_k )
%
% and the Accelerated Gradient Descent (AGD) 
%
% u_{k+1} = y_k - alpha*Ep( y_k )
%
% y_{k+1} = u_{k+1} + ((lambda-1)/(lambda+1))*( u_{k+1} - u_k );
%
% methods for a stongly convex and (locally) Lipschitz smooth 
% objective E with derivative Ep.
%
% Input
%   E : the objective function
%   Ep : the derivative of the objective function
%   mu : (an estimate of) the strong convexity constant for E
%   L : (an estimate of) the (local) Lipschitz smoothness 
%       constant for E
%   MaxIts : the maximal number of iterations
%   uinit : an initial guess for the minimizer
%
% Output
%   energyPSD(1:MaxIts+1) : the value of the objective at every 
%                           iteration of PSD
%   energyAGD(1:MaxIts+1) : the value of the objective at every 
%                           iteration of AGD
%   upsd : the approximate minimizer after MaxIts iterations 
%          of PSD
%   uagd : the approximate minimizer after MaxIts iterations 
%          of AGD
%

% We first run PSD
  alphaPSD = 0.5/L;

  energyPSD = zeros(MaxIts+1,1);
  energyPSD(1) = E( uinit );

  upsd = uinit;
  for i=1:MaxIts
    upsd = PSD( alphaPSD, upsd, Ep );
    energyPSD(i+1) = E( upsd );
  end

% Next we run AGD
  alphaAGD = 0.5/L;
  lambdaAGD = sqrt( 2*L/mu );

  energyAGD = energyPSD;
  uagd = uinit;
  yagd = uinit;
  for i=1:MaxIts
    [uagd, yagd ] = AGD(alphaAGD, lambdaAGD, uagd, yagd, Ep );
    energyAGD(i+1) = E( uagd );
  end
end


function u = PSD( alpha, uold, Ep )
% This function does one iteration of the PSD scheme
%
% Input
%   alpha : the line search parameter
%   uold : the previous approximate minimizer
%   Ep : the derivative of the objective
%
% Output
%   u : the next approximate minimizer
%
  u = uold - alpha.*Ep( uold );
end

function [u, y] = AGD(alpha, lambda, uold, yold, Ep )
% This function does one iteration of the AGD scheme
%
% Input
%   alpha : the line search parameter
%   lambda : the extrapolation parameter
%   uold : the previous approximate minimizer
%   yold : the previous extrapolated value
%   Ep : the derivative of the objective
%
% Output
%   u : the next approximate minimizer
%   y : the next approximate extrapolation
%
  u = yold - alpha.*Ep( yold );
  y = u + ((lambda-1.)/(lambda+1.)).*( u - uold );
end
