function c = array2cell(a)
  [~,~,H] = size(a);
  c = {};
  for i = 1:H
    c{i} = a(:,:,i);
  end
