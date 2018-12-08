slice1=19;
slice2=25;
%Load primary eigenvectors:
Vx=fib.dir0(1,:,:,slice1:slice2);
Vy=fib.dir0(2,:,:,slice1:slice2);
Vz=fib.dir0(3,:,:,slice1:slice2);

Vx=reshape(Vx,96,96,7);
Vy=reshape(Vy,96,96,7);
Vz=reshape(Vz,96,96,7);


% f1=size(fib.dir0,1);
% f2=size(fib.dir0,2);
% f3=size(fib.dir0,3);
f1=src.dimension(1);
f2=src.dimension(2);
f3=src.dimension(3);
[X,Y,Z]=meshgrid(1:f2,1:f1,slice1:slice2);        

%Choose a starting point from the dataset
x1=trk.tracts(1,1);
y1=trk.tracts(2,1);
z1=trk.tracts(3,1);
x2=trk.tracts( 1,trk.length(1));
y2=trk.tracts( 2,trk.length(1));
z2=trk.tracts( 3,trk.length(1));
    
x_path=[x1,x2]; 
y_path=[y1,y2];
z_path=[z1,z2];
i1=round(x_path(1));
j1=round(y_path(1));
z1=round(z_path(1));
i2=round(x_path(2));
j2=round(y_path(2));
z2=round(z_path(2));

z1=6; %6 because this track starts from slice 24 which is the 6th volume for Dijkstra
z2=6;
%Dijkstra
[d,pi,pj,pz,QS]=dijkstra_function_3D_version2(i1,j1,z1,X,Y,Z,Vy,Vx,Vz,fib.fa0(:,:,slice1:slice2));

%Backtracking and Plotting
figure;h1=plot3(X(i1,j1,z1),Y(i1,j1,z1),Z(i1,j1,z1),'o','MarkerSize',10); hold on;
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);
      
xroma=cmap(mod(1,size(cmap,1))+1,:);                           
[path_lin,path_col,path_vol,B_backtrack]=dijkstra_backtrack_fun_3D(d,pi,pj,pz,i2,j2,z2); 
idx=sub2ind(size(X),path_lin,path_col,path_vol);    
plot3(X(i2,j2,z2),Y(i2,j2,z2),Z(i2,j2,z2),'p','MarkerSize',10,'Color',xroma);
h=plot3(X(idx),Y(idx),Z(idx),'-','Linewidth',2, 'Color',xroma);
drawnow;

%clear slice1 slice2 Vx Vy Vz f1 f2 f3 x1 y1 z1 x2 y2 z2 x_path y_path z_path i1 j1 z1 i2 j2 z2 xroma X Y Z cmap 