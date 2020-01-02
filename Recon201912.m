close all;
clear all;
clc

addpath('./UniversalFunctions');
addpath('./Read&Write');
addpath('./jRL_Deconvolution');
addpath('./Read&Write');

folderpath = uigetdir('G:\Data\');
disp(folderpath);
disp('Loading config and experiment data.');
config = ReadConfig(folderpath);

tstr = datestr(now, 30);
filepath = fullfile(config.Reconstruction.ResultPath,tstr);
if ~exist(filepath,'dir')
  mkdir(filepath);
end

Files = dir(fullfile(config.Data.folderpath,'*.tif'));
FileNames = {Files.name};

filename = fullfile(config.Data.folderpath, FileNames{1});
Info = imfinfo(filename);
Slice = size(Info,1);


OTFConfig = otf_gen(config);
m = size(OTFConfig,1);
OTF = zeros(m*2, m*2);
OTF(m/2:m*3/2-1,m/2:m*3/2-1) = OTFConfig;
imwrite(OTF, fullfile(filepath, 'OTF.bmp'));
[M,N] = size(OTF);
RLimit = get_otf_paraments(OTF);
x = linspace(1,M,M);
cent = fix(M/2+1);
[X,Y] = meshgrid(x,x);
Rv = abs((X-cent) + 1j*(Y-cent));

alpha = 0.1*RLimit;
beta = 0.5;
gamma = alpha;
delta = 1;
OSpara = 0.03*RLimit;
OSR = 20;

MaskF = Rv < RLimit;
MaskH = Rv < 0.5*RLimit;
MaskD = Rv < 2*RLimit;
Mask2 = exp(-Rv.^2/(2*gamma*gamma));
%Mask2 = Rv < 0.2*RLimit;

Crop = config.Data.isCrop;
locx = config.Data.CropStart(1);
locy = config.Data.CropStart(2);
c_size = config.Data.CropSize;

WFCell = {};
OSCell = {};
SRCell = {};
SROSCell = {};
OSSRCell = {};
WFOSCell = {};
WFDeCell = {};

z = 1;
for z = 1:fix(Slice/4)
    DataCell = {};

    i = 1 + 4*(z-1);
    for j = 1:4
    Img = im2double(imread(filename,i+j-1));
    if Crop
      DataCell{j} = Img(locx:locx+c_size-1, locy:locy+c_size-1);
    else
      DataCell{j} = Img;
    end
    DataCell{j} = DoubleSize(DataCell{j});
    end

    WriteTIFF(DataCell, fullfile(filepath, ['RawData_z',int2str(z),'.tif']));

    Data = zeros(M,N,4);
    Is = zeros(1,4);
    for j = 1:4
    Data(:,:,j) = DataCell{j};
    %temp = abs(ifft2(MaskF.*fftshift(fft2(temp))));
    %Data(:,:,j) = temp / max(temp(:));
    Is(j) = sum(sum(DataCell{j}));
    end
    Is = Is / max(Is(1));
    %imwrite(Data(:,:,1)/max(max(Data(:,:,1))), fullfile(filepath, ['WF_z',int2str(z),'.bmp']));


    OSImg = OC1(Data(:,:,2:4), OTF, OSpara, OSR);

    %figure, imshow(OSImg);
    %WriteTIFF(array2cell(OSImg),fullfile(filepath,'OSImg.tif'));
    OSDeImg = jRL_deconvolve(OSImg, ones(M,N), OTF, 50);
    %WriteTIFF(array2cell(OSDeImg),fullfile(filepath,'OSDeImg.tif'));
    

    WFOS = CombineImg2(Data(:,:,1), OSImg, alpha, 0.8);
    %figure,imshow(WFOS);
    %WriteTIFF(array2cell(WFOS),fullfile(filepath,'WFOS.tif'));

    WFDe = jRL_deconvolve(Data(:,:,1), ones(M,N), OTF, 50);
    %figure,imshow(WFDe);
    %WriteTIFF(array2cell(WFDe), fullfile(filepath, 'WFDe.tif'));

    %[SIPattern, ~, ~, ~] = Estimate(Data(:,:,2:4), OTF, Data(:,:,1));
    %WriteTIFF(array2cell(SIPattern),fullfile(filepath, 'SIPattern.tif'));

    fm = zeros(size(Data));
    for n = 1:4
        f = fftshift(fft2(Data(:,:,n)));
        fm(:,:,n) = 20*log(abs(f));
    end
    %WriteTIFF(array2cell(fm),fullfile(filepath,'FM.tif'));



    [Phase, K, MfList] = IlluParam(Data(:,:,2:4),OTF,Data(:,:,1));

%{
    r = RLimit;
    xc = 0;
    yc = xc;

    theta = linspace(0,2*pi);
    x = r*cos(theta) + xc;
    y = r*sin(theta) + yc;

    figure, hold on,
    plot(x,y),
    scatter(K(1,2),K(1,1),'r'),
    scatter(K(2,2),K(2,1),'g'),
    scatter(K(3,2),K(3,1),'b'),
    axis equal
    hold off;
%}
    
    MF = [0, 0.05, 0.2, 0.3];
    MF(2:4) = MfList(2,:);
    %I0 = [1, 0.5, 1, 0.8];
    I0 = Is;
    SIPattern = ones(M,N,3);
    for n = 1:3
        SIPattern(:,:,n) = get_pattern2(K(n,:),Phase(n),M, MF(n+1), I0(n+1));
    end
    Pattern = ones(M,N,4)*I0(1);
    Pattern(:,:,2:4) = SIPattern;


    %WriteTIFF(array2cell(Pattern),fullfile(filepath,'Pattern.tif'));

    SRImg = jRL_deconvolve(Data, Pattern, OTF, 50);
    %WriteTIFF(array2cell(SRImg), fullfile(filepath, 'SRImg.tif'));
    
    SROS = CombineImg2(SRImg, OSImg, alpha, beta);
    SROSDe = CombineImg2(SRImg, OSDeImg, alpha, beta);
    %WriteTIFF(array2cell(SROS), fullfile(filepath, 'SROS.tif'));
    %WriteTIFF(array2cell(SROSDe), fullfile(filepath, 'SROSDe.tif'));

    OF = Data(:,:,1) - delta * OSImg;
    OF(OF<0) = 0;
    OF = abs(ifft2(fftshift(fft2(OF)).*Mask2));
    
    DataIF = zeros(size(Data));
    DataIF(:,:,1) = SROS;
    
    for j = 2:4
        OStemp = SSOS(Data(:,:,j), OTF, OSpara, OSR);
        OFtemp = Data(:,:,1)*Is(j) - delta * OStemp;
        OFtemp(OFtemp < 0) = 0;
        DataIF(:,:,j) = Data(:,:,j) - abs(ifft2(fftshift(fft2(OFtemp)).*Mask2));
    end
    DataIF(DataIF < 0) = 0;

    WriteTIFF(array2cell(DataIF), fullfile(filepath,['DataIF_z',int2str(z), '.tif']))
    
    [PhaseIF, KIF, MfListIF] = IlluParam(DataIF(:,:,2:4),OTF,DataIF(:,:,1));
    MF = [0, 0.05, 0.2, 0.3];
    MF(2:4) = MfListIF(2,:);
    %I0 = [1, 0.5, 1, 0.8];
    I0 = Is;
    SIPatternIF = ones(M,N,3);
    for n = 1:3
        SIPatternIF(:,:,n) = get_pattern2(KIF(n,:),PhaseIF(n),M, MF(n+1), I0(n+1));
    end
    PatternIF = ones(M,N,4)*I0(1);
    PatternIF(:,:,2:4) = SIPatternIF;
    tic
    SRIFImg = jRL_deconvolve(DataIF, PatternIF, OTF, 50);
    %WriteTIFF(array2cell(SRIFImg), fullfile(filepath, 'SRIFImg.tif'));
    toc
    SSSR = CombineImg2(SRIFImg, OSImg, alpha, beta);
    %WriteTIFF(array2cell(SSSR), fullfile(filepath, 'SSSR.tif'));
    
    WFCell{z} = Data(:,:,1);
    OSCell{z} = OSImg;
    OSDeCell{z} = OSDeImg;
    SRCell{z} = SRImg;
    SROSCell{z} = SROS;
    SROSDeCell{z} = SROSDe;
    OSSRCell{z} = SRIFImg;
    WFOSCell{z} = WFOS;
    WFDeCell{z} = WFDe;
    
    disp(z);
    
    WriteTIFF(WFCell, fullfile(filepath, 'zWFCell.tif'));
    WriteTIFF(OSCell, fullfile(filepath, 'zOSCell.tif'));
    WriteTIFF(OSDeCell, fullfile(filepath, 'zOSDeCell.tif'));
    WriteTIFF(WFOSCell, fullfile(filepath, 'zWFOSCell.tif'));
    WriteTIFF(SRCell, fullfile(filepath, 'zSRCell.tif'));
    WriteTIFF(SROSCell, fullfile(filepath, 'zSROSCell.tif'));
    WriteTIFF(SROSDeCell, fullfile(filepath, 'zSROSDeCell.tif'));
    WriteTIFF(OSSRCell, fullfile(filepath, 'zOSSRCell.tif'));
    WriteTIFF(WFDeCell, fullfile(filepath, 'zWFDeCell.tif'));
    
end

    WriteTIFF(WFCell, fullfile(filepath, 'zWFCell.tif'));
    WriteTIFF(OSCell, fullfile(filepath, 'zOSCell.tif'));
    WriteTIFF(OSDeCell, fullfile(filepath, 'zOSDeCell.tif'));
    WriteTIFF(WFOSCell, fullfile(filepath, 'zWFOSCell.tif'));
    WriteTIFF(SRCell, fullfile(filepath, 'zSRCell.tif'));
    WriteTIFF(SROSCell, fullfile(filepath, 'zSROSCell.tif'));
    WriteTIFF(SROSDeCell, fullfile(filepath, 'zSROSDeCell.tif'));
    WriteTIFF(OSSRCell, fullfile(filepath, 'zOSSRCell.tif'));
    WriteTIFF(WFDeCell, fullfile(filepath, 'zWFDeCell.tif'));
%{
temp = cell2array(WFCell);
temp = temp / max(temp(:))*255;
WFCell2 = array2cell(temp);
WriteTIFF(WFCell2, fullfile(filepath, 'zWFCell2.tif'));


temp = cell2array(OSCell);
temp = temp / max(temp(:))*255;
OSCell2 = array2cell(temp);
WriteTIFF(OSCell2, fullfile(filepath, 'zOSCell2.tif'));

temp = cell2array(WFOSCell);
temp = temp / max(temp(:))*255;
WFOSCell2 = array2cell(temp);
WriteTIFF(WFOSCell2, fullfile(filepath, 'zWFOSCell2.tif'));

temp = cell2array(SRCell);
temp = temp / max(temp(:))*255;
SRCell2 = array2cell(temp);
WriteTIFF(SRCell2, fullfile(filepath, 'zSRCell2.tif'));

temp = cell2array(SROSCell);
temp = temp / max(temp(:))*255;
SROSCell2 = array2cell(temp);
WriteTIFF(SROSCell2, fullfile(filepath, 'zSROSCell2.tif'));

temp = cell2array(OSSRCell);
temp = temp / max(temp(:))*255;
OSSRCell2 = array2cell(temp);
WriteTIFF(OSSRCell2, fullfile(filepath, 'zOSSRCell2.tif'));

temp = cell2array(WFDeCell);
temp = temp / max(temp(:))*255;
WFDeCell2 = array2cell(temp);
WriteTIFF(WFDeCell2, fullfile(filepath, 'zWFDeCell2.tif'));
    %}