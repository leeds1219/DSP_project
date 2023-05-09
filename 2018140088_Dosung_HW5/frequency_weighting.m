function xx = frequency_weighting(fc,f,v)

xx = exp(-((log2(f)-log2(fc)).^2)/2*(v^2));

end