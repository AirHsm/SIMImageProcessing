close all
clear all
clc



m = 512;
z = 512;

obj = zeros(m, m, z);
r = 9;

cent = fix([m/2+1,m/2+1,z/2+1]);

a = 5;
b = 5;
c = 64;

x = linspace(1,m,m)-cent(1);
y = linspace(1,m,m)-cent(2);
z = linspace(1,z,z)-cent(3);

[X,Y,Z] = meshgrid(x,y,z);
Rv = sqrt(X.^2+Y.^2+Z.^2);

step = 30;

for t = 0:10:180
  theta = deg2rad(t);
  %theta = pi/2;
  for h = 0:10:360
    phi = deg2rad(h);
    u = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    temp = (X.*u(1)+Y.*u(2)+Z.*u(3))./norm(u,2)./Rv;
    temp = rad2deg(acos(temp));
    obj = (abs(temp)<(2)) + obj;
  end
end
%{
theta = pi/2;
phi = 0;
u = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
temp = (X.*u(1)+Y.*u(2)+Z.*u(3))./norm(u,2)./Rv;
temp = rad2deg(acos(temp));
obj = (abs(temp)<(5));
%}
obj(obj>1) = 1;
obj = obj.*(Rv<200);
figure,imshow(obj(:,:,cent(3)));
s = 3;
save('obj_512x512.mat','obj');
% img = imgaussfilt3(obj,s, 'FilterDomain', 'frequency');
%{
h3 = fspecial3('gaussian',9,s);
h2 = fspecial('gaussian',9,s);
data3 = imfilter(obj,h3,'same');
data2 = imfilter(obj,h2,'same');
figure,
subplot(221), montage(obj(:,:,cent(3)-4:cent(3)+4));
subplot(222), montage(data2(:,:,cent(3)-4:cent(3)+4));
subplot(223), montage(data3(:,:,cent(3)-4:cent(3)+4));
%subplot(224), montage(h3(:,:,cent(3)-4:cent(3)+4));

figure,
subplot(221), montage(obj);
subplot(222), montage(data2);
subplot(223), montage(data3);
subplot(224), montage(h3);
%}