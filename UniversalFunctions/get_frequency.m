function k_e = get_frequency(sim_img, otf)
  f = abs(fftshift(fft2(sim_img)));

  m = size(sim_img,1);
  cent = fix(m/2+1);
  x = linspace(1,m,m);
  y = linspace(1,m,m);
  [X,Y] = meshgrid(x,y);

  r = sqrt((X-cent).^2+(Y-cent).^2);
  R_limit = get_otf_paraments(otf);

  Zone = (r > 0.6*R_limit).*(r < 0.8*R_limit).* ( Y <= cent) .*(r <= R_limit).*(abs(Y-cent)>1).*(abs(X-cent)>1);
  temp = f.*Zone;

  %k = zeros(1,2);
  [~, Index] = max(temp(:));
  [k_row,k_col] = ind2sub(size(temp),Index);
  k0 = [cent,cent]-[k_row,k_col];
  %starting optimize

  %FreqOpt0 = @(k)FreqOpt(k,sim_img);
  %options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
  %k_e = fminunc(FreqOpt0, k0, options);

  %options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',1000,'MaxIter',500,'Display','notify');
  % [k_e,~] = fminsearch(FreqOpt0,k,options);

  k_e = k0;
end
