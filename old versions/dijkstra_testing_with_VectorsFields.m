% clear all
close all
if 1>0
%Efarmogi tou djk se VectorFields
%Epilogi ton simeion me xrisi tou getpts


% [x_test,y_test]=meshgrid(-10:0.1:10,-10:0.1:10);
% figure(3);quiver(x_test,y_test,sin(x_test),sin(2*y_test));

%Dimiourgia ton Vector fields

[X,Y]=meshgrid(-2:0.02:2,-2:0.02:2);
Vx=-Y;
Vy=X;

Vx=sin(X);
Vy=sin(2*Y);
[Vx,Vy]=gradient(peaks(X,Y));
Cl=sqrt(Vx.^2+Vy.^2);
Cl=Cl/max(Cl(:));
else
    load 'real_data_slice15';
    Vx=x;
    Vy=y;
    [X,Y]=meshgrid(1:256,1:256);
    Cl=cl_volume15;
end

figure(1);quiver(X,Y,Vx,Vy);
% figure(2);quiver(x,y,sin(x),sin(2*y));
%figure(2);quiver(x,y,sin(x),cos(2*y));
fig=figure(1);

hold on;
%[x_path,y_path]=getpts(fig); 

 x_path=[30,10];
 y_path=[20,20];
% x_t=110;
% y_t=120;
i1=round(x_path(1));
j1=round(y_path(1));
i2=round(x_path(end));
j2=round(y_path(end));



%metatropi ton dianismaton apo cartesianes se polikes :
%angle_th_eigvctr einai i gonia se rad
%radius_eigvctr to mikos ton dianismaton
[angle_th_eigvctr,radius_eigvctr]=cart2pol(Vx,Vy);
axis equal

[d,pi,pj,QS]=dijkstra_testing_with_VectorsFields_djk_fun(i1,j1,X,Y,Vx,Vy,Cl);

h1=plot(X(i1,j1),Y(i1,j1),'o','MarkerSize',10);
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);
dests=randint(100,2,[1+2,size(X,1)-2]);

for k=1:100
    xroma=cmap(mod(k,size(cmap,1))+1,:);
    i2=dests(k,1);
    j2=dests(k,2);
    [path_lin,path_col,b]=dijkstra_testing_with_VectorsFields_backtrack_fun(d,pi,pj,i2,j2);
    idx=sub2ind(size(X),path_lin,path_col);    
    plot(X(i2,j2),Y(i2,j2),'p','MarkerSize',5,'Color',xroma);
    h=plot(X(idx),Y(idx),'-','Linewidth',2, 'Color',xroma); %get gca get(h)
    drawnow;
end

%quiver(y0,x0,A1(index_cl),A2(index_cl));
%plot(pathy,pathx,'r-','Linewidth',2);
%plot(pathy,pathx,'r-');