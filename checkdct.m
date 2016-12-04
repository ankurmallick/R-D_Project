function y=checkdct(x)
L=size(x,3);
y=zeros(size(x));
x0=0;
for i=1:L
    x0=x0+x(:,:,i);
    y(:,:,i)=idct2(x0);
end
end