function pattern = get_pattern(k,phase,m,mf,I0)
  % pattern = get_pattern(k,phase,m,mf,I0)
  % modified at 20181119.
  % I0 is the mean of the raw image where 0<I0<1.
  cent = fix(m/2+1);
  x = linspace(1,m,m);
  y = linspace(1,m,m);
  [X,Y] = meshgrid(x,y);

  pattern = (1+mf*cos(2*pi*(k(2)/m.*(X-cent)+k(1)/m.*(Y-cent))+phase))*I0;
  pattern = pattern/max(pattern(:))*I0;
end
