%%% Dimiourgia synthetikon data gia elegxo tou 3D dijkstra algorithm

%% dimiourgia tis 3D kampilis gia to ground truth
[X,Y,Z]=meshgrid(-2:0.0627:2,-2:0.0627:2,-2:0.0627:2);
fig1=figure; %plot(X(:),Y(:),'.');
axis equal;
% [xp,yp]=getpts;
xp=[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
yp=[-1.8000,-0.9442,-0.5354,-0.2239,0.2239,0.9832,1.4000];

zp=linspace(-1.5,1.5,length(xp));
N=(xp(end)-xp(1))/0.1+1;    
xx=linspace(xp(1),xp(end),N);
yy=spline(xp,yp,xx);
zz=spline(xp,zp,xx);

%% fiber1 
%kataskeui idiodianismaton:
Vx=zeros(size(X));
Vy=zeros(size(X));
Vz=zeros(size(X));
B=zeros(size(X));
for k=1:length(xx)-1
    %vriskoume to pio kontino stn kampili pou ftiaksame parapano
    z=( X(:)-xx(k) ).^2 + ( Y(:)-yy(k) ).^2 + ( Z(:)-zz(k) ).^2; 
    [minz,mink]=min(z);
    
    B(mink)=1;
    Vx(mink)=xx(k+1)-xx(k);
    Vy(mink)=yy(k+1)-yy(k);
    Vz(mink)=zz(k+1)-zz(k);
end

hold on;
plot3(xp,yp,zp,'-or');
plot3(xx,yy,zz,'-ok');
%quiver3(X,Y,Z,Vx,Vy,Vz);

%% kataskeuazoume tin maska tou inpainting algorithm(paper) diffusion
mask=ones(3,3,3)/26;
mask(2,2,2)=0;

%% fiber 2
% % % xp2=-[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
% % % yp2=yp;
% % % N=(xp2(end)-xp2(1))/0.1+1;
% % % xx2=linspace(xp2(1),xp2(end),abs(N ) )
% % % yy2=spline(xp2,yp2,xx2);
% % % 
% % % Vx2=zeros(size(X));
% % % Vy2=zeros(size(X));
% % % for k=1:length(xx2)-1
% % %     z=(X(:)-xx2(k)).^2 + (Y(:)-yy2(k)).^2;
% % %     [minz,mink]=min(z);
% % %     Vx2(mink)=xx2(k+1)-xx2(k);
% % %     Vy2(mink)=yy2(k+1)-yy2(k);
% % % end
% % %  hold on;
% % %  plot(xp2,yp2,'-or');
% % %  plot(xx2,yy2,'-ok');
% % % quiver(X,Y,Vx2,Vy2);
% % % Cl2=sqrt(Vx2.^2+Vy2.^2);
% % % Cl2=Cl2/max(Cl2(:));

%% Diffusion me tin maska tou inpainting algorithm(paper) stis 3 diastaseis

for iter=1:40
    Vz=convn(Vz,mask,'same');
    Vx=convn(Vx,mask,'same');
    Vy=convn(Vy,mask,'same');
    
end

%quiver3(X,Y,Z,Vx,Vy,Vz);
Cl=sqrt(Vx.^2+Vy.^2+Vz.^2);
Cl=Cl/max(Cl(:));
%h=quiver(X,Y,Vx,Vy,'m');
%hold on;
%quiver(X,Y,Vx2,Vy2,'m');

%% NOISE:
% Vx_noise=Vx(:)+0.0027;
% Vy_noise=Vy(:)+0.0027 ;
% Vx_noise=reshape(Vx_noise,41,41);
% Vy_noise=reshape(Vy_noise,41,41); 
% figure;quiver(X,Y,Vx_noise,Vy_noise,'m'); title('with noise');

