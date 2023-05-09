% 4. A Musical Illusion
fs = 22050;

C = zeros(1,length(key2note(40)));
for i = 1:9
    C = key2note(40+8*(5-i)) + C;
end

C_sharp = zeros(1,length(key2note(41)));
for i = 1:9
    C_sharp = key2note(41+8*(5-i)) + C_sharp;
end

D = zeros(1,length(key2note(42)));
for i = 1:9
    D = key2note(42+8*(5-i)) + D;
end

D_sharp = zeros(1,length(key2note(43)));
for i = 1:9
    D_sharp = key2note(43+8*(5-i)) + D_sharp;
end

E = zeros(1,length(key2note(44)));
for i = 1:9
    E = key2note(44+8*(5-i)) + E;
end

F = zeros(1,length(key2note(45)));
for i = 1:9
    F = key2note(45+8*(5-i)) + F;
end

F_sharp = zeros(1,length(key2note(46)));
for i = 1:9
    F_sharp = key2note(46+8*(5-i)) + F_sharp;
end

G = zeros(1,length(key2note(47)));
for i = 1:9
    G = key2note(47+8*(5-i)) + G;
end

G_sharp = zeros(1,length(key2note(48)));
for i = 1:9
    G_sharp = key2note(48+8*(5-i)) + G_sharp;
end

A = zeros(1,length(key2note(49)));
for i = 1:9
    A = key2note(49+8*(5-i)) + A;
end

A_sharp = zeros(1,length(key2note(50)));
for i = 1:9
    A_sharp = key2note(50+8*(5-i)) + A_sharp;
end

B = zeros(1,length(key2note(51)));
for i = 1:9
    B = key2note(51+8*(5-i)) + B;
end

xx = cat(2,C,C_sharp);
xx = cat(2,xx,D);
xx = cat(2,xx,D_sharp);
xx = cat(2,xx,E);
xx = cat(2,xx,F);
xx = cat(2,xx,F_sharp);
xx = cat(2,xx,G);
xx = cat(2,xx,G_sharp);
xx = cat(2,xx,A);
xx = cat(2,xx,A_sharp);
xx = cat(2,xx,B);
xx = repmat(xx,1,5);

f = figure;
spectrogram(xx,fs);
f2 = figure;
plotspec(xx,fs);

soundsc(xx,fs)