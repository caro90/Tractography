%%% Dimiourgia synthetikon data gia elegxo tou 2D dijkstra algorithm

% %% dimiourgia tis 2D kampilis gia to ground truth
% %[X,Y]=meshgrid(1:100,1:100);
% [X,Y]=meshgrid(-2:0.0628:2,-2:0.0628:2);
% %figure; %plot(X(:),Y(:),'.');
% axis equal;
% % [xp,yp]=getpts;
% xp=[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
% yp=[-1.8000,-0.9442,-0.5354,-0.2239,0.2239,0.9832,1.4000];
% N=(xp(end)-xp(1))/0.1+1;    
% xx=linspace(xp(1),xp(end),N);
% yy=spline(xp,yp,xx);
% 
% %% fiber1 
% %kataskeui idiodianismaton:
% Vx=zeros(size(X));
% Vy=zeros(size(X));
% 
% %Ston B kratame to track tis gnostis kabilis 
% %gia na to sigkrinoume me to track tou dijkstra (ground truth )
% B=zeros(size(X));
% 
% for k=1:length(xx)-1
%     %vriskoume to pio kontino simeio stin kampili
%     %z=abs( X(:)-xx(k) ) + abs(Y(:)-yy(k)); 
%     z=( X(:)-xx(k) ).^2 + (Y(:)-yy(k)).^2; 
%     [minz,mink]=min(z);
%     B(mink)=1;
%     Vx(mink)=xx(k+1)-xx(k);
%     Vy(mink)=yy(k+1)-yy(k);
% end
% 
% %% kataskeuazoume tin maska tou inpainting algorithm(paper) diffusion
% mask=ones(3,3)/8;
% mask(2,2)=0;
% 
% %% Diffusion fiber 1
% 
% for iter=1:50
%     Vx=conv2(Vx,mask,'same');
%     Vy=conv2(Vy,mask,'same');
% end
% 
% %% Plotting fiber 1
% % figure;hold on;
% % plot(xp,yp,'-or');
% % plot(xx,yy,'-ok');
% % quiver(X,Y,Vx,Vy);title('fiber1');
% % Cl=sqrt(Vx.^2+Vy.^2);
% % Cl=Cl/max(Cl(:));
%  
%  
% %% fiber 2
% xp2=-[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
% yp2=yp;
% N=(xp2(end)-xp2(1))/0.1+1;
% xx2=linspace(xp2(1),xp2(end),abs(N ) );
% yy2=spline(xp2,yp2,xx2);
% 
% Vx2=zeros(size(X));
% Vy2=zeros(size(X));
% for k=1:length(xx2)-1
%       z=(X(:)-xx2(k)).^2 + (Y(:)-yy2(k)).^2;
%       %z=abs(X(:)-xx2(k)) + abs(Y(:)-yy2(k));
%       [minz,mink]=min(z);
%       B(mink)=1;
% %     Vx2(mink)=xx2(k)-xx2(k+1);
% %     Vy2(mink)=yy2(k)-yy2(k+1);
%   Vx2(mink)=xx2(k+1)-xx2(k);
%   Vy2(mink)=yy2(k+1)-yy2(k);
% end
% 
% %% Diffusion fiber 2
% for iter=1:50
%     Vx2=conv2(Vx2,mask,'same');
%     Vy2=conv2(Vy2,mask,'same');
% end
% 
% %% Plotting fiber 2
% % figure;hold on;
% % plot(xp2,yp2,'-or');
% % plot(xx2,yy2,'-ok');
% % quiver(X,Y,Vx2,Vy2);title('fiber 2');
% % Cl2=sqrt(Vx2.^2+Vy2.^2);
% % Cl2=Cl2/max(Cl2(:));
% 
% 
% 
% %% prosthesi ton idiodianismaton tou Fiber 1 kai Fiber 2.Ypologismos sinolikou Cl
% Vx=Vx+Vx2;
% Vy=Vy+Vy2;
% Cl=sqrt(Vx.^2+Vy.^2);
% Cl=Cl/max(Cl(:));
% 
% 
% %% plotting fiber1 + fiber 2
% figure;hold on;
% plot(xp,yp,'-or');
% plot(xx,yy,'-ok');
% h=quiver(X,Y,Vx,Vy);
% plot(xp2,yp2,'-or');
% plot(xx2,yy2,'-ok');
% 
% %% NOISE testing:
% % Vx_noise=Vx(:)+0.0027;
% % Vy_noise=Vy(:)+0.0027 ;
% % Vx_noise=reshape(Vx_noise,41,41);
% % Vy_noise=reshape(Vy_noise,41,41);
% % 
% % figure;quiver(X,Y,Vx_noise,Vy_noise,'m'); title('with noise');
% 
% %% other types of vector fields
% % [X,Y]=meshgrid(-2:0.05:2,-2:0.05:2);
% %  %Vx=cos(X^2+Y^2);
% %  %Vy=X-Y^2+1;
% %  Vx1=-Y;
% %  Vy1=X;
% %  Vx2=exp(Y);
% %  Vy2=X;
% % %Vx=X.^3;
% % %Vy=X.^3;
% % Vx=Vx1+Vx2;
% % Vy=Vy1+Vy2;
% % 
% % Cl=sqrt(Vx.^2+Vy.^2);
% % Cl=Cl/max(Cl(:));
% % figure;quiver(X,Y,Vx,Vy);
% % hold on;
% 
% 
% %% Ground truth checking
% % 
% % D = bwdist(B);
% % x=mean( D(:).*B_backtrack(:));
% % g_t=x/ans;
% 



%% new fiber1
N=128;
[X,Y]=meshgrid(1:N,1:N);
Vax=rand(128)./10;
Vay=rand(128)./10;
Vbx=rand(128)./10;
Vby=rand(128)./10;

Vay(:,N/2)=1;
Vbx(N/2,:)=1;

Vax(:,N/2)=0;
Vby(N/2,:)=0;
%%
Ima=zeros(N);
Ima(:,N/2)=1;
dista=bwdist(Ima);

Imb=zeros(N);
Imb(N/2,:)=1;
distb=bwdist(Imb);

La=( N/2 -dista )/ (N/2);
Lb=( N/2 -distb )/ (N/2);

%% kataskeuazoume tin maska tou inpainting algorithm(paper) diffusion
mask=ones(3,3)/8;
mask(2,2)=0;

%% Diffusion fiber 1

for iter=1:10
    Vax=conv2(Vax,mask,'same');
    Vay=conv2(Vay,mask,'same');
    Vbx=conv2(Vbx,mask,'same');
    Vby=conv2(Vby,mask,'same');
end


% figure;quiver(Vax,Vay);
% hold on;
% quiver(Vbx,Vby,'color',[1 0 0]);

% V1x=max(Vax,Vbx);
% V1y=max(Vay,Vby);
% V1z=zeros(N);
% 
% V2x=min(Vax,Vbx);
% V2y=min(Vay,Vby);
% V2z=zeros(N);

L1=zeros(N);
L2=zeros(N);
L3=zeros(N);

% La=sqrt(Vax.^2+Vay.^2);
% Lb=sqrt(Vbx.^2+Vby.^2);

idx=find(La>Lb);
V1x(idx)=Vax(idx);
V1y(idx)=Vay(idx);
V2x(idx)=Vbx(idx);
V2y(idx)=Vby(idx);
L1(idx)=La(idx);
L2(idx)=Lb(idx);

idx=find(La<=Lb);
V1x(idx)=Vbx(idx);
V1y(idx)=Vby(idx);
V2x(idx)=Vax(idx);
V2y(idx)=Vay(idx);
L1(idx)=Lb(idx);
L2(idx)=La(idx);

V1x=reshape(V1x,N,N);
V1y=reshape(V1y,N,N);
V2x=reshape(V2x,N,N);
V2y=reshape(V1y,N,N);

%L3=zeros(1,N*N);

Cl=(L1-L2)./L1;
Cp=(L2-L3)./L1;
Cs=L3./L1;

figure;quiver(V1x,V1y);
hold on;
quiver(V2x,V2y,'color',[1 0 0]);

% figure;imshow(Cl,[]);title('cl');
% figure;imshow(Cp,[]);title('cp');
% 
% FA=( (L1-L2).^2+(L1-L3).^2+(L2-L3).^2 )./(L1.^2 + L2.^2 +L3.^2);
% figure;imshow(FA,[]);title('Fa');



