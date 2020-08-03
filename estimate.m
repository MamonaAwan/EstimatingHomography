%reading image
Img1=imread('img1.jpg');
Img2=imread('img2.jpg');
I1 =im2double(rgb2gray(Img1));
I2 =im2double(rgb2gray(Img2));
[Row,Col]=size(I1);
[Corners1x Corners1y]=CornersDetect(I1);
[Corners2x Corners2y]=CornersDetect(I2);
r1=size(Corners1x,1);
r2=size(Corners2x,1);

Match1x=zeros(450,1);
Match1y=zeros(450,1);
WindowSize=11;
Win=(WindowSize-1)/2;
for i=1+Win:1:r2-Win
    MinNCC=0.0;
    Li=0; Lj=0;
    for a=1+Win:1:r1-Win
        CurrentNCC=0.0;
        NumNCC=0.0;     DenNCC=0.0;
        RdenNCC=0.0;    LdenNCC=0.0;
        Lmean=0.0;      Rmean=0.0;
        for x=-Win:1:Win
            for y=-Win:1:Win
                Lmean=Lmean+I1(Corners1x(a)+x,Corners1y(a)+y);
                Rmean=Rmean+I2(Corners2x(i)+x,Corners2y(i)+y);
            end
        end
        Lmean=Lmean/(WindowSize^2);
        Rmean=Rmean/(WindowSize^2);
        for x=-Win:1:Win
            for y=-Win:1:Win
                NumNCC=NumNCC+((I1(Corners1x(a)+x,Corners1y(a)+y)-Lmean)*(I2(Corners2x(i)+x,Corners2y(i)+y)-Rmean));
                RdenNCC=RdenNCC+((I2(Corners2x(i)+x,Corners2y(i)+y)-Rmean)^2);
                LdenNCC=LdenNCC+((I1(Corners1x(a)+x,Corners1y(a)+y)-Lmean)^2);
            end
        end
        DenNCC=sqrt(RdenNCC*LdenNCC);
        CurrentNCC=NumNCC/DenNCC;
        if MinNCC<CurrentNCC
            MinNCC=CurrentNCC;
            Li=Corners1x(a);
            Lj=Corners1y(a);
        end
    end
    Match1x(i)=Li;
    Match1y(i)=Lj;
end

c1=[Match1y Match1x];
c2=[Corners2y Corners2x];
sz=size(c1,1);

%% Ransac Algoritm
npts=12;
minerror=1e10;
m = fix((1-.1)*sz)-20;
for j=1:1:450
    i=randi(sz,npts,1);
    x=c1(i,1); y=c1(i,2);
    u=c2(i,1); v=c2(i,2);
    % Taking ones and zeros to add desired columns
    o= ones(size(x));
    z= zeros(size(x));
 
    %Constructing A Matrix and Calculating all parts
    A_part1= [x y o z z z -u.*x -u.*y -u];
    A_part2= [z z z x y o -v.*x -v.*y -v];
    A=[A_part1; A_part2];
    [~,~,V] = svd(A,0);
    h = V(:,end);
    H = reshape(h,3,3)';
    d2 = zeros(sz,1); 
    for k = 1:sz
        xt = H*[c1(k,:),1]'; 
        t = xt/xt(3);
        xt = t(1:2);
        d = c2(k,:)-xt';
        d2(k) = sqrt(d(1)^2+d(2)^2);
    end
    [Y I] = sort(d2);
    if sum(Y(1:sz)) < minerror
        minerror = sum(Y(1:sz));
        inliers = I(1:m);
        outliers = I(m+1:end);
    end 
end
p=4;

Im=[Img1 Img2];
figure;
imshow(Im);
hold on;
for i = inliers'
    plot([c1(i,1),Col+c2(i,1)],[c1(i,2),c2(i,2)],'-og','linewidth',1);
end
for i = outliers'
    plot([c1(i,1),Col+c2(i,1)],[c1(i,2),c2(i,2)],'-or','linewidth',1);
end

x=c1(inliers,1); y=c1(inliers,2);
u=c2(inliers,1); v=c2(inliers,2);
x1=x(1:p); y1=y(1:p);
u1=x(1:p); v1=x(1:p);
o= ones(size(x1));
z= zeros(size(x1));
A_part1= [x1 y1 o z z z -u1.*x1 -u1.*y1 -u1];
A_part2= [z z z x1 y1 o -v1.*x1 -v1.*y1 -v1];
A=[A_part1; A_part2];
[~,~,V] = svd(A);
h1 = V(:,9);
H1 = reshape(h1,3,3)';
H=RefineHomo(c1,c2,H1);
Hinv=inv(H);
Hinv=Hinv./Hinv(3,3);

%% method 1 (only one method is used at a time, other is commented)
% using maketransform
T=maketform('projective',[u(1:4) v(1:4)],[x(1:4) y(1:4)]);
[im2t,xdataim2t,ydataim2t]=imtransform(Img2,T);
%xdataim2t and ydataim2t store the bounds of the transformed img2
xdataout=[min(1,xdataim2t(1)) max(size(Img1,2),xdataim2t(2))];
ydataout=[min(1,ydataim2t(1)) max(size(Img1,1),ydataim2t(2))];
%transform both images with the computed xdata and ydata
im2t=imtransform(Img2,T,'XData',xdataout,'YData',ydataout);
im1t=imtransform(Img1,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);
ims=[im1t/2+im2t/2];
imd=uint8(abs(double(im1t)-double(im2t)));
ims=max(im1t,im2t);
figure,imshow(ims);

%% method 2 (only one method is used at a time, other is commented)
% using bilinear interplation
[M N C]=size(Img1);
h = H'; h = h(:);
c14 = Ft(h,[1,1,N,1,1,M,N,M]);
e = [1,3,5,7];
f = [2,4,6,8];
xmin = round(min(c14(e)));xmax = round(max(c14(e)));
ymin = round(min(c14(f)));ymax = round(max(c14(f)));
img = zeros(875,1125,C);
Img3=[zeros(750,625,3) Img2(:,350:1000,:)];
img = alignImg(img,Img3,H,xmin+15,ymin-20);
img = alignImg(img,Img1(:,1:570,:),eye(3,3),xmin,ymin);
figure;
imshow(img);
