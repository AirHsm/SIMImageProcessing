function coef =  FreqOpt(k, sim_img)
  m = size(sim_img,1);
  cent = fix(m/2+1);
  x = linspace(1,m,m);
  y = linspace(1,m,m);
  [X,Y] = meshgrid(x,y);

  S = exp(-1j*2*pi.*(k(2)/m*(X-cent)+k(1)/m*(Y-cent)))*sim_img;
  Sf = fft2(S);
  coef = -abs(sum(sum(fft2(sim_img).*conj(Sf)))/(sum(sum(Sf.*conj(Sf)))+eps));
end
