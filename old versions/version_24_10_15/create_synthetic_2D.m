%%% Dimiourgia synthetikon data gia elegxo tou 2D dijkstra algorithm

%% dimiourgia tis 2D kampilis gia to ground truth
[X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
figure; %plot(X(:),Y(:),'.');
axis equal;
% [xp,yp]=getpts;
xp=[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
yp=[-1.8000,-0.9442,-0.5354,-0.2239,0.2239,0.9832,1.4000];
N=(xp(end)-xp(1))/0.1+1;    
xx=linspace(xp(1),xp(end),N);
yy=spline(xp,yp,xx);

%% fiber1 
%kataskeui idiodianismaton:
Vx=zeros(size(X));
Vy=zeros(size(X));

%Ston B kratame to track tis gnostis kabilis 
%gia na to sigkrinoume me to track tou dijkstra
B=zeros(size(X));

for k=1:length(xx)-1
    %vriskoume to pio kontino simeio stin kampili
    z=abs( X(:)-xx(k) ) + abs(Y(:)-yy(k)); 
    %z=( X(:)-xx(k) ).^2 + (Y(:)-yy(k)).^2; 
    [minz,mink]=min(z);
    B(mink)=1;
    Vx(mink)=xx(k+1)-xx(k);
    Vy(mink)=yy(k+1)-yy(k);
end

%% kataskeuazoume tin maska tou inpainting algorithm(paper) diffusion
mask=ones(3,3)/8;
mask(2,2)=0;

%% Diffusion fiber 1

for iter=1:20
    Vx=conv2(Vx,mask,'same');
    Vy=conv2(Vy,mask,'same');
end

%% Plotting fiber 1
hold on;
plot(xp,yp,'-or');
plot(xx,yy,'-ok');
%quiver(X,Y,Vx,Vy);
Cl=sqrt(Vx.^2+Vy.^2);
Cl=Cl/max(Cl(:));
 
 
%% fiber 2
xp2=-[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
yp2=yp;
N=(xp2(end)-xp2(1))/0.1+1;
xx2=linspace(xp2(1),xp2(end),abs(N ) );
yy2=spline(xp2,yp2,xx2);

Vx2=zeros(size(X));
Vy2=zeros(size(X));
for k=1:length(xx2)-1
    z=(X(:)-xx2(k)).^2 + (Y(:)-yy2(k)).^2;
    [minz,mink]=min(z);
    Vx2(mink)=xx2(k+1)-xx2(k);
    Vy2(mink)=yy2(k+1)-yy2(k);
end

%% Diffusion fiber 2

for iter=1:20
    Vx2=conv2(Vx2,mask,'same');
    Vy2=conv2(Vy2,mask,'same');
end

%% Plotting fiber 2
hold on;
plot(xp2,yp2,'-or');
plot(xx2,yy2,'-ok');
% quiver(X,Y,Vx2,Vy2);
Cl2=sqrt(Vx2.^2+Vy2.^2);
Cl2=Cl2/max(Cl2(:));


%% prosthesi ton idiodianismaton tou Fiber 1 kai Fiber 2.Ypologismos sinolikou Cl
Vx=Vx+Vx2;
Vy=Vy+Vy2;
Cl=sqrt(Vx.^2+Vy.^2);
Cl=Cl/max(Cl(:));

%% plotting fiber1 + fiber 2
h=quiver(X,Y,Vx,Vy);



%% NOISE testing:
% Vx_noise=Vx(:)+0.0027;
% Vy_noise=Vy(:)+0.0027 ;
% Vx_noise=reshape(Vx_noise,41,41);
% Vy_noise=reshape(Vy_noise,41,41);
% 
% figure;quiver(X,Y,Vx_noise,Vy_noise,'m'); title('with noise');

