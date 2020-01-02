clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is very important file.
% reconstruction for simulation in the paper 


% dependence
addpath '../UniversalFunctions'
addpath '../jRL_Deconvolution'
addpath '../Read&Write'

filepath = 'E:\workspace\projects\SIMImaging\MatlabSIM\3DImagingSimulation\ExperimentData\20191008';

load TestData.mat
H = 9;
Data = zeros(M,M,H);
Pattern = zeros(size(Data));
OTF = abs(fftshift(fft2(psf(:,:,1))));

RLimit = get_otf_paraments(OTF);
x = linspace(1,M,M);
cent = fix(M/2+1);
[X,Y] = meshgrid(x,x);
Rv = abs((X-cent) + 1j*(Y-cent));

WF = BluredData{1};
WriteTIFF(array2cell(WF),fullfile(filepath,'WF.tif'));

WFDe = jRL_deconvolve(WF, ones(size(WF)), OTF, 50);
WriteTIFF(array2cell(WFDe), fullfile(filepath,'WFDe.tif'));


for i = 1:H
    Data(:,:,i) = BluredData{i+1};
    Pattern(:,:,i) = Illu{i+1};
end

WriteTIFF(array2cell(Data),fullfile(filepath,'Data.tif'));
WriteTIFF(array2cell(Pattern),fullfile(filepath,'Pattern.tif'));

OSImg = zeros(size(WF));
for i = 1:3:9
    OSImg = OSImg + get_OS(Data(:,:,i),WF)/3;
end

figure,imshow(OSImg/ max(OSImg(:)));
WriteTIFF(array2cell(OSImg), fullfile(filepath,'OSImg.tif'));
OSDe = jRL_deconvolve(OSImg, ones(size(OSImg)), OTF, 50);
figure,imshow(OSDe/max(OSDe(:)));
WriteTIFF(array2cell(OSDe), fullfile(filepath,'OSDe.tif'));
%{
Data4 = zeros(M,M,4);
Data4(:,:,1) = WF;
Data4(:,:,2:4) = Data(:,:,1:3:9);
SI4 = zeros(size(Data4));
SI4(:,:,1) = Illu{1};
SI4(:,:,2:4) = Pattern(:,:,1:3:9);

SRjRL4 = jRL_deconvolve(Data4, SI4, OTF, 50);
figure,imshow(SRjRL4/max(SRjRL4(:)));
%}
SRjRL = jRL_deconvolve(Data, Pattern, OTF, 50);
figure,imshow(SRjRL/max(SRjRL(:)));
WriteTIFF(array2cell(SRjRL),fullfile(filepath,'SRjRL.tif'));

%SROS = CombineImg2(SRjRL, OSDe, 0.3*RLimit,0.5);
%figure,imshow(SROS/max(SROS(:)));

%a = 0.3;
b = 0.5;
for a = 0:0.1:1
    SROS = CombineImg2(SRjRL, OSDe, a*RLimit,b);
    SROS = (SROS-min(SROS(:))/(max(SROS(:))-min(SROS(:))));
    WriteTIFF(array2cell(SROS),fullfile(filepath,['SROS_a_',num2str(a),'_b_',num2str(b),'.tif']));
end

a = 0.3;
for b = 0:0.1:1
    SROS = CombineImg2(SRjRL, OSDe, a*RLimit,b);
    %SROS = (SROS-min(SROS(:))/(max(SROS(:))-min(SROS(:))));
    WriteTIFF(array2cell(SROS),fullfile(filepath,['SROS_a_',num2str(a),'_b_',num2str(b),'.tif']));
end



delta = 1;
for g = 0:0.1:1
    DataIF = zeros(size(Data));
    gamma = g*RLimit;
    %Mask2 = exp(-Rv.^2/(2*gamma*gamma));
    Mask2 = Rv < gamma;
    for i = 1:9
        OStemp = get_OS(Data(:,:,i), WF);
        OFtemp = Data(:,:,i) - OStemp;
        OFtemp(OFtemp < 0) = 0;
        temp = abs(ifft2(fftshift(fft2(OFtemp)).*Mask2));
        DataIF(:,:,i) = Data(:,:,i) - delta*temp;
    end
    OSSR = jRL_deconvolve(DataIF, Pattern, OTF, 50);
    WriteTIFF(array2cell(OSSR), fullfile(filepath,['OSSR_g_',num2str(g),'_d_',num2str(delta),'.tif']));
end


g = 0.3;
gamma = g*RLimit;
Mask2 = Rv < gamma;
for delta = 0:0.1:1
    DataIF = zeros(size(Data));
    for i = 1:9
        OStemp = get_OS(Data(:,:,i), WF);
        OFtemp = Data(:,:,i) - OStemp;
        OFtemp(OFtemp < 0) = 0;
        temp = abs(ifft2(fftshift(fft2(OFtemp)).*Mask2));
        DataIF(:,:,i) = Data(:,:,i) - delta*temp;
    end
    OSSR = jRL_deconvolve(DataIF, Pattern, OTF, 50);
    WriteTIFF(array2cell(OSSR), fullfile(filepath,['OSSR_g_',num2str(g),'_d_',num2str(delta),'.tif']));
end
%figure,imshow(OSSR/max(OSSR(:)));
%{
Data4IF = zeros(M,M,4);
Data4IF(:,:,1) = OSImg;
Data4IF(:,:,2:4) = DataIF(:,:,1:3:9);

OSSR4 = jRL_deconvolve(Data4IF, SI4, OTF, 50);
figure,imshow(OSSR4/max(OSSR4(:)));
%}
%}


