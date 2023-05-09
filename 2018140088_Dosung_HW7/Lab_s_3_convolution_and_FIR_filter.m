
% 2.1
% dconvdemo

% 2.2
% a
load echart.mat
f1 = figure; imshow(echart)
% b
yyl = firfilt(echart(65,:),[1 -1]);
f2=figure;
subplot(2,1,1)
stem(echart(65,:))
title('original row')

subplot(2,1,2)
stem(yyl)
title('filtered row')
% c

% d

% Additional problem
load CT_chest.mat

yy_0 = firfilt(chest(1,:),[1 -1]);
for m = 1:442
    yy = firfilt(chest(m,:),[1 -1]);
    yy_0 = cat(1,yy_0,yy);
end

xx_0 = firfilt(chest(:,1),[1 -1]);
for n = 1:440
    xx = firfilt(chest(:,n),[1 -1]);
    xx_0 = cat(1,xx_0,xx);
end

f3 = figure;
subplot(2,2,1)
imshow(chest)
title('original')

subplot(2,2,2)
imshow(transpose(xx_0))
title('x edge')

subplot(2,2,3)
imshow(yy_0)
title('y edge')

subplot(2,2,4)
imshow(transpose(xx_0)+yy_0)
title('x and y edge')