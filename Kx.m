function caps = Kx(Kmax, Kwidth, xmax, patchnum, traits)

store = zeros(patchnum, length(traits));
for k = 1:patchnum
    store(k,:) = Kmax(1,k)*exp((-(traits-xmax(1,k)).^2/(2*(Kwidth(1,k)).^2)));
end

caps = store;