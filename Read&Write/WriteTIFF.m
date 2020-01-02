%Write tiff file from an array MxNxH
function WriteTIFF(data, filename)
  if exist(filename)
    delete(filename)
  end
  Slice = size(data,2);
  for i = 1:Slice
    temp = data{i};
    temp = temp / max(temp(:));
    img = uint8(temp*256);
    imwrite(img,filename,'WriteMode','Append');
  end
