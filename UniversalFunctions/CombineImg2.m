function Result = CombineImg2(SR, OS,a,b)
  [M,~] = size(SR);
  x = linspace(1,M,M);
  y = linspace(1,M,M);
  [x,y] =meshgrid(x,y);
  cent = M/2+1;
  Rv = abs((x-cent)+1i*(y-cent));

  %Mask = exp(-Rv.^2/(2*a*a));
  Mask = Rv < a;

  Result = abs(ifft2(fftshift((1-b)*(1-Mask).*fftshift(fft2(SR))+b*Mask.*fftshift(fft2(OS)))));
end
