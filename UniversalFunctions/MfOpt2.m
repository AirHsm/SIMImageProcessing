% MfOpt2
function res = MfOpt2(a, wf, k, phi, si)
  [M,N] = size(wf);

  x = linspace(1,N,N);
  y = linspace(1,M,M);
  [X,Y] = meshgrid(x,y);
  cent = M/2+1;

  % temp = wf.*(a(1)+a(2)*cos(2*pi*(k(2)/M.*(X-cent)+k(1)/M.*(Y-cent))+phi));
  temp = wf.*(1+a*cos(2*pi*(k(2)/M.*(X-cent)+k(1)/M.*(Y-cent))+phi));
  res = sum(sum(conj(temp-si).*(temp-si)));
