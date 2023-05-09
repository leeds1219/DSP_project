function xx = key2note(X, keynum, dur)
    j = sqrt(-1);
    f = 440*2^((keynum-49)/12); 
    fs = 11025;
    tt = 0:(1/fs):dur;
    xx = real(X.*exp(j*2*pi*f*tt));
end