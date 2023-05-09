function xx = key2note(keynum)
    j = sqrt(-1);
    f = 440*2^((keynum-49)/12); 
    fs = 22050; % freq sample
    dur = 0.5; % duration of each note
    tt = 0:(1/fs):dur;
    v = 2; % the variance
    fc = 440;
    X = exp(-((log2(f)-log2(fc)).^2)/2*(v^2));
    xx = real(X.*exp(j*2*pi*f*tt));
end