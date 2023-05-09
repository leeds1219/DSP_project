% Demo for aliasing
% 1 run con2dis.m

% Compare Frequency of the component
% 2
tt = 0:0.001:1;
a1 = signal(11,tt);
a2 = signal(7,tt);
a3 = signal(2,tt);
% 3 
b1 = 10.*a1+3.*a2+a3; % high frequency components?
b2 = a1+3.*a2+10.*a3;
f1 = figure;
subplot(2,1,1)
plot(tt,b1)
subplot(2,1,2)
plot(tt,b2)
% b1 이 high frequency component 가 강하다
% Restoration of signal
load unknown_signal.mat

tt = 0:1/100:0.1; % fs 100 Hz

Fsfun = fscreate(tt,unkn3,10,'zoh');

aa = Fsfun('coef');
bb = Fsfun('period');

f2 = figure;
subplot(2,1,1)
stem(abs(aa))
subplot(2,1,2)
plot(tt,unkn3)

bb = Fsfun('phase');

temp = sort(abs(aa),'descend');
temp2 = unique(temp);

% 이걸 1000 Hz fs 로 복원 do not use interpolation

n = 0:1/1000:100;
restored_signal = 0.63662*2*cos(((2*pi*10/1000)*n)+pi) + 0.491816*2*cos((2*pi*50/1000)*n);

f3 = figure;
plot(n,restored_signal)

% Prepare a test image
% 4
xx = imread('abdomen.bmp');
% 5
xx_max = max(max(xx)); % 255
xx_min = min(min(xx)); % 0
% 6
f4 = figure;imshow(xx)
% 7
% What parts in the image mainly contain high frequency components 
% and low frequency components respectively? 
% Think about objects in the edges of each
% high contrast가 반복되는 부분이 high frequency components 일 것이다

% Down-sample the signal
% 8 hint: xx(1:interval:end,1:interval:end)
interval = 3;
xx_down = xx(1:interval:end,1:interval:end);

f5 =figure; imshow(xx_down)

% the high frequency component part ha aliasing effect
% 흑백의 대조가 많은 부분에서 생략되어 

% Reconstruct the original siganl from the down-sampled signal
% 9  zoh, hint : interp2(……. nearestz)
[rows,cols,~] = size(xx_down);
[X,Y]=meshgrid(1:cols,1:rows);
[Xq,Yq] = meshgrid(1:1/interval:cols,1:1/interval:rows);
recontructed_image = interp2(X,Y,double(xx_down),Xq,Yq,'nearest');

recontructed_image = uint8(recontructed_image);% convert back to uint8

f6 = figure; imshow(recontructed_image);

% 10 linear interpolation, hint: interp2(…..hlineari)

[rows, cols, ~] = size(xx_down);
[X, Y] = meshgrid(1:cols, 1:rows);
[Xq, Yq] = meshgrid(1:1/interval:cols, 1:1/interval:rows);
reconstructed_image_linear = interp2(X, Y, double(xx_down), Xq, Yq, 'linear');
reconstructed_image_linear = uint8(reconstructed_image_linear);

f7 = figure; imshow(reconstructed_image_linear);
% 11 zoh 와 linear interpolation 비교 설명 

