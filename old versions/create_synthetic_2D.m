clear all;
% close all;
[X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
figure; plot(X(:),Y(:),'.');
axis equal;
% [xp,yp]=getpts;
xp=[-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000];
yp=[-1.8000,-0.9442,-0.5354,-0.2239,0.2239,0.9832,1.4000];
N=(xp(end)-xp(1))/0.1+1;
xx=linspace(xp(1),xp(end),N)
yy=spline(xp,yp,xx);
Vx=zeros(size(X));
Vy=zeros(size(X));
for k=1:length(xx)-1
    z=(X(:)-xx(k)).^2 + (Y(:)-yy(k)).^2;
    [minz,mink]=min(z);
    Vx(mink)=xx(k+1)-xx(k);
    Vy(mink)=yy(k+1)-yy(k);
end

hold on;
plot(xp,yp,'-or');
plot(xx,yy,'-ok');
quiver(X,Y,Vx,Vy);

mask=ones(3,3)/8;
mask(2,2)=0;
Cl=sqrt(Vx.^2+Vy.^2);
Cl=Cl/max(Cl(:));

%% fiber 2
xp2=-fliplr([-1.6000,-1.1973,-0.2044,0.8274,1.3336,1.6938,1.8000]);
yp2=yp;
N=(xp2(end)-xp2(1))/0.1+1;
xx2=linspace(xp2(1),xp2(end),N)
yy2=spline(xp2,yp2,xx2);
Vx2=zeros(size(X));
Vy2=zeros(size(X));
for k=1:length(xx2)-1
    z=(X(:)-xx2(k)).^2 + (Y(:)-yy2(k)).^2;
    [minz,mink]=min(z);
    Vx2(mink)=xx2(k+1)-xx2(k);
    Vy2(mink)=yy2(k+1)-yy2(k);
end
%%

% hold on;
% plot(xp2,yp2,'-or');
% plot(xx2,yy2,'-ok');
% quiver(X,Y,Vx2,Vy2);
% mask=ones(3,3)/8;
% mask(2,2)=0;
% 
% Cl=sqrt(Vx2.^2+Vy2.^2);
% Cl=Cl/max(Cl(:));

for iter=1:50
    Vx=conv2(Vx,mask,'same');
    Vy=conv2(Vy,mask,'same');
end
quiver(X,Y,Vx,Vy,'m');
