%Author: Ankur Mallick
%Framewise tomographic video reconstruction
close all
clear all
clc

%For Video Data
load('frames.mat');
%For CT slice data
%load(‘ctframes.mat'); %Uncomment
% N=144;
% Nframes=1;
% I=imresize(imread('lena.jpg'),[N,N]);

x_best=zeros(N,N,Nframes);
x_init=x_best;
x_den=x_init;
rad=[];
y=[];
rad_init=[];
m1=min(I(:));
m2=max(I(:));
for i=1:Nframes
    thetas=rem(i,3):3:179;
    % Defining function handles
    A=@(x)(dct2proj(x,thetas));
    AT=@(x)(dct2projT(x,thetas,N));
    D=@(x)(checkdct(x));
    DI=@(x)(checkdctI(x));
    DT=@(x)(checkdctT(x));
    DTI=@(x)(checkdctTI(x));
    % Video corruption
    x_best(:,:,i)=10*mat2gray(I(:,:,i),[m1,m2]); %Low intensity image (0 to 10)
    d_best=DI(x_best(:,:,i));
    rad(:,:,i)=A(d_best);
    if(sum(sum(rad(:,:,i)<0))>0)
        rad(:,:,i)=rad(:,:,i)-min(min(rad(:,:,i)));
    end
    y=poissrnd(rad(:,:,i)); %Poisson corrupted projections
    L=10; %Patch length
    b=10^-10;
    lambda=0.001;
    % Initialisation
    rad_init(:,:,i)=global_denoise(y,L); %Globally denoised projections
    x0=iradon(rad_init(:,:,i),thetas,N); %FBP is used here for the initial guess but not in the actual algorithm
    if(sum(sum(x0<0))>0)
        %Some negative elements
        x0=x0-min(x0(:));
    end
    x_init(:,:,i)=x0;
    d_init=DI(x_init(:,:,i));
    d=d_init;
    alpha=1; %stepsize
    endflag=0;
    Ad=A(d);
    T=Ad-y.*log(Ad+b);
    S=sum(sum(abs(d)));
    F=sum(T(:))+lambda*S(:);
    F_init=F;
    disp('hi');
    %% Denoising
    for j=1:2000
        dir=AT(1-y./(Ad+b)); %Gradient
        while(1)
            %Minimising differentiable part
            d1=d-alpha*dir;
            d2=wthresh(d1,'s',lambda*alpha);
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
    iter=j;
    disp('Number of iterations: ');
    disp(iter);
    disp('Value of step size(exit condition): ');
    disp(alpha);
    x_den(:,:,i)=D(d);
    R(:,:,i)=Ad;
    disp(i);
end
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
    %     title(strcat('Stage 1 denoised frame ',num2str(i)));
    figure
    imshow(mat2gray(x_den(:,:,i),[m,M]));
    title(strcat('Stage 2 denoised frame ',num2str(i)));
    %     figure
    %     imshow(mat2gray(x_best(:,:,i),[m,M]));
    %     title(strcat('Original frame ',num2str(i)));
end
E1=norm(x_init(:)-x_best(:))/norm(x_best(:));
disp('RMSE of initial guess: ');
disp(E1);
E2=norm(x_den(:)-x_best(:))/norm(x_best(:));
disp(‘Final RMSE: ');
disp(E2);save('data_vidframe.mat','x_den','x_best','x_init','rad','rad_init','R');