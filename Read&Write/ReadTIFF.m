%Read tiff file into an array (MxNxH)
function output = ReadTIFF(filename)
  Info = imfinfo(filename);
  tif = 'tif';
  format = Info.Format;
  if strcmp(format, tif) == 0
    disp('This is not a .tif file, please check.')
  end
  Slice = size(Info,1);
  Width = Info.Width;
  Height = Info.Height;

  output = zeros(Height, Width, Slice);

  for i = 1 : Slice
    output(:,:,i) = im2double(imread(filename, i));
  end
