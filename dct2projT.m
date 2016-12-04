function y=dct2projT(x,thetas,N)
L=size(thetas,1);
y0=dct2(iradon(x(:,:,1),thetas(1,:),'linear','none',1,N));
S=[size(y0), L];
y=zeros(S);
y(:,:,1)=y0;
if(L>1)
    for i=2:L
        y0=dct2(iradon(x(:,:,i),thetas(i,:),'linear','none',1,N));
        y(:,:,1:i)=y(:,:,1:i)+y0(:,:,ones(1,i));
    end
end
end