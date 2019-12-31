function a = cell2array(c)
  H = max(size(c));
  [M,N] = size(c{1});
  a = zeros(M,N,H);

  for i = 1:H
    a(:,:,i) = c{i};
  end
