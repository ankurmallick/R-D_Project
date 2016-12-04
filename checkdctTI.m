function y=checkdctTI(x)
%Upper triangular matrix
L=size(x,3);
y=zeros(size(x));
x0=0;
for i=L:-1:1
    y(:,:,i)=idct2(x(:,:,i)-x0);
    x0=x(:,:,i);
end
end