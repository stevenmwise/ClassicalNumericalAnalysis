function y = WilkinsonFact(x)
  y = 1.;
  for i=1:20
      y = y .*( x .- i );  
  end
end
