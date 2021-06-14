% The method of Cs to determine the order of consistency of a linear multistep method.
% Input
% a : The coefficients of the first characteristic polynomial
% b : The coefficients of the second characteristic polynomial
% m : The number of C that one wants to compute
% Output
% res: The number C_m. A method is consistent to order exactly p if C_0 = ... = C_p = 0, but C_{p+1} != 0

function res = methodCs( a, b, m )
  if m == 0
    res = sum( a );
  else
    res = 0.;
    q = length(a);
    factmminusone = factorial(m-1);
    factm = m*factmminusone;
    for j=1:q
      res = res + a(j)*(j-1)^m/factm - b(j)*(j-1)^(m-1)/factmminusone;
    end
  end
end
