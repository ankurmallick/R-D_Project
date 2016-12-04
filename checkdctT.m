function y=checkdctT(x)
L=size(x,3);
y=zeros(size(x));
for i=1:L
    y0=dct2(x(:,:,i));
    y(:,:,1:i)=y(:,:,1:i)+y0(:,:,ones(1,i));
end
end