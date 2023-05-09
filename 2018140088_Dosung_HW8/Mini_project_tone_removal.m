% Mini project tone removal
[xx,fs] = audioread('SunshineSquare.wav');
xx = xx';
f1=figure;
specgram(xx, [], fs) 

% nulling filters
% 1.
w1 = 2*pi*4715/fs;
A1 = -2*cos(w1);
hh1 = [1,A1,1];
y1= filter(hh1,1,xx);

ww = -pi:pi/100:pi;
H1 = freqz(hh1,1,ww);
f2=figure;
subplot(2,1,1)
specgram(y1,[],fs)
subplot(2,1,2)
plot(ww,abs(H1)); 
% 2.
w2 = 2*pi*3150/fs;
A2 = -2*cos(w2);
hh2 = [1,A2,1];
y2 = filter(hh2,1,xx);

ww = -pi:pi/100:pi;
H2 = freqz(hh2,1,ww);
f3=figure;
subplot(2,1,1)
specgram(y2,[],fs)
subplot(2,1,2)
plot(ww,abs(H2)); 
% 3.
w3 = 2*pi*1570/fs;
A3 = -2*cos(w3);
hh3 = [1,A3,1];
y3 = filter(hh3,1,xx);

ww = -pi:pi/100:pi;
H3 = freqz(hh3,1,ww);
f4=figure;
subplot(2,1,1)
specgram(y3,[],fs)
subplot(2,1,2)
plot(ww,abs(H3)); 
%4.
w4 = 0;
A4 = -2*cos(w4);
hh4 = [1,A4,1];
y4 = filter(hh4,1,xx);

ww = -pi:pi/100:pi;
H4 = freqz(hh4,1,ww);
f5=figure;
subplot(2,1,1)
specgram(y4,[],fs)
subplot(2,1,2)
plot(ww,abs(H4));
%5
w5 = 2*pi*4720/fs;
A5 = -2*cos(w5);
hh5 = [1,A5,1];

%6
w6 = 2*pi*1600/fs;
A6 = -2*cos(w6);
hh6 = [1,A6,1];

% 1-6 은 함수를 하나 만들어서 input argument 에 w를 받아주면 
% 코드를 더 간편하게 정리할  수 있다.

hh = conv(hh1,hh2);
hh = conv(hh,hh3);
hh = conv(hh,hh4);
hh = conv(hh,hh5);
hh = conv(hh,hh6);
yy = filter(hh,1,xx);
f0=figure; specgram(yy,[],fs)
player = audioplayer(yy,fs);
play(player)