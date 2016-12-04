function y=dct2proj(x,thetas)
L=size(thetas,1);
xprev=x(:,:,1);
y0=radon(idct2(xprev),thetas(1,:));
S=[size(y0), L];
y=zeros(S);
y(:,:,1)=y0;
if(L>1)
    for i=2:L
        xprev=xprev+x(:,:,i);
        y(:,:,i)=radon(idct2(xprev),thetas(i,:));
    end
end
end