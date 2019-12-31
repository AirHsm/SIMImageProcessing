function output = DoubleSize(img)
  [m,n] = size(img);
  f = fftshift(fft2(img));
  f_d = zeros(m*2,n*2);
  s_x = fix(m/2+1);
  s_y = fix(m/2+1);
  f_d(s_x:s_x+m-1,s_y:s_y+n-1) = f;
  output = abs(ifft2(f_d));
