function y=checkdctI(x)
%Lower triangular matrix
L=size(x,3);
y=zeros(size(x));
x0=0;
for i=1:L
    y(:,:,i)=dct2(x(:,:,i)-x0);
    x0=x(:,:,i);
end
end