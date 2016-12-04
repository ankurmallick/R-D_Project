%Author: Ankur Mallick
%Joint tomographic video reconstruction
close all
clear all
clc

load('frames.mat');
% N=144;
% Nframes=1;
% I=imresize(imread('lena.jpg'),[N,N]);

%% Defining function handles
%N=256;
thetas=randi([0,179],[Nframes,60]); %projections for frames
%Random angles do not work well with multiple (5) frames
for i=1:Nframes
    j=rem(i,3);
    thetas(i,:)=j:3:179;
end
A=@(x)(dct2proj(x,thetas));
AT=@(x)(dct2projT(x,thetas,N));
D=@(x)(checkdct(x));
DI=@(x)(checkdctI(x));
DT=@(x)(checkdctT(x));
DTI=@(x)(checkdctTI(x));

%% Video corruption
x_best=zeros(size(I));
m1=min(I(:));
m2=max(I(:));
for i=1:Nframes
x_best(:,:,i)=10*mat2gray(I(:,:,i),[m1,m2]); %Low intensity image (0 to 10)
end
% for i=1:Nframes
%     R=radon(x_best(:,:,i),thetas(i,:));
%     disp(min(min(R)));
% end
d_best=DI(x_best);
rad=A(d_best);
rad(rad<0)=0;
% disp(min(min(rad)));
% if(sum(rad(:)<0)>0)
%     rad=rad-min(rad(:));
% end
y=poissrnd(rad); %Poisson corrupted projections
L=10; %Patch length

%% Initialisation
rad_init=y;
x_init=zeros(N,N,Nframes);
for i=1:Nframes
    rad_init(:,:,i)=global_denoise(y(:,:,i),L); %Globally denoised projections
    x0=iradon(rad_init(:,:,i),thetas(i,:),N); %FBP is used here for the initial guess but not in the actual algorithm
    if(sum(sum(x0<0))>0)
        %Some negative elements
        x0=x0-min(x0(:));
    end
    %x0(x0<0)=0;
    x_init(:,:,i)=x0; 
end
d_init=DI(x_init);
d=d_init;
lambda=zeros(1,Nframes);
lambda(1)=0.001; %regularizer
lambda(2:Nframes)=0.1*lambda(1);
%Later on try with different lambdas
alpha=1; %stepsize
endflag=0;
b=10^-10;
Ad=A(d);
T=Ad-y.*log(Ad+b);
S=sum(sum(abs(d)));
F=sum(T(:))+lambda*S(:);
F_init=F;
T_best=rad-y.*log(rad+b);
S_best=sum(sum(abs(d_best)));
F_best=sum(T_best(:))+lambda*S_best(:);

%% Denoising
for i=1:2000
    disp(F);
    dir=AT(1-y./(Ad+b)); %Gradient
    while(1)
        %Minimising differentiable part
        d1=d-alpha*dir;
        d2=zeros(size(d1));
        for j=1:Nframes
            %ISTA
            d2(:,:,j)=wthresh(d1(:,:,j),'s',lambda(j)*alpha);
        end
        %Imposing non-negativity constraint
        d3=-DTI(d2);
        d3(d3<0)=0;
        d_new=d2+DT(d3);
        Ad_new=A(d_new);
        T_new=Ad_new-y.*log(Ad_new+b);
        S_new=sum(sum(abs(d_new)));
        F_new=sum(T_new(:))+lambda*S_new(:);
        if(F_new<F)
            %Updating values
            F=F_new;
            d=d_new;
            Ad=Ad_new;
            alpha=1.5*alpha;
            break;
        end
        if(alpha<10^-30)
            %Exit condition
            endflag=1;
            break;
        end
        alpha=0.75*alpha;
    end
    if(endflag==1)
        break;
    end
end
iter=i;
disp('Number of iterations: ');
disp(iter);
disp('Value of step size(exit condition): ');
disp(alpha);
x_den=D(d);
R=Ad;
MSE1=sum(((x_init(:)-x_best(:)).^2))/numel(x_best);
PSNR1=20*log10(max(x_best(:)))-10*log10(MSE1);
disp('PSNR of initial guess: ')
disp(PSNR1);
MSE2=sum(((x_den(:)-x_best(:)).^2))/numel(x_best);
PSNR2=20*log10(max(x_best(:)))-10*log10(MSE2);
disp(‘Final PSNR: ')
disp(PSNR2);
m=min([min(x_best(:)),min(x_den(:)),min(x_init(:))]);
M=max([max(x_best(:)),max(x_den(:)),max(x_init(:))]);
for i=1:Nframes
%     figure
%     imshow(mat2gray(x_init(:,:,i),[m,M]));
%     title(strcat(‘Initial denoised frame ',num2str(i)));
    figure
    imshow(mat2gray(x_den(:,:,i),[m,M]));
    title(strcat(‘Final denoised frame ',num2str(i)));
%     figure
%     imshow(mat2gray(x_best(:,:,i),[m,M]));
%     title(strcat('Original frame ',num2str(i)));
end
E1=norm(x_init(:)-x_best(:))/norm(x_best(:));
disp('RMSE of initial guess: ');
disp(E1);
E2=norm(x_den(:)-x_best(:))/norm(x_best(:));
disp(‘Final RMSE: ');
disp(E2);
save('data_vid.mat','x_den','x_best','x_init','rad','rad_init','R');