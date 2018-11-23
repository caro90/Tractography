%Load primary eigenvectors:
slice1=19;
slice2=25;
Vx=fib.dir0(:,:,slice1:slice2,1);
Vy=fib.dir0(:,:,slice1:slice2,2);
Vz=fib.dir0(:,:,slice1:slice2,3);

f1=size(fib.dir0,1);
f2=size(fib.dir0,2);
f3=size(fib.dir0,3);
%f1 f2 reverse
[X,Y,Z]=meshgrid(1:f2,1:f1,slice1:slice2);        

%Choosing a starting point from the dataset
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

%z=1 because we restrict our slices
z1=6;
z2=6;
[d,pi,pj,pz,QS]=dijkstra_function_3D_version2(i1,j1,z1,X,Y,Z,Vx,Vy,Vz,fib.fa0(:,:,19:25));

%Backtracking and Plotting
h1=plot3(X(i1,j1,z1),Y(i1,j1,z1),Z(i1,j1,z1),'o','MarkerSize',10); hold on;
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);
for k=1:1      
    xroma=cmap(mod(k+1,size(cmap,1))+1,:);                           
    [path_lin,path_col,path_vol,B_backtrack]=dijkstra_backtrack_fun_3D(d,pi,pj,pz,i2,j2,z2); 
    idx=sub2ind(size(X),path_lin,path_col,path_vol);    
    plot3(X(i2,j2,z2),Y(i2,j2,z2),Z(i2,j2,z2),'p','MarkerSize',10,'Color',xroma);
    figure;
    h=plot3(X(idx),Y(idx),Z(idx),'-','Linewidth',2, 'Color',xroma);
    drawnow;

    %plotting test gia real data
%     if (flag==3)
%             figure;imshow(Cl(:,:,1),[0 1]);hold on;
%             for i=1:length(path_vol) 
% 
%                 if eq(path_vol(i),1)
%                     idx1=idx(i);
%                     plot3(X(idx1),Y(idx1),Z(idx1),'-','Linewidth',2, 'Color',xroma);
%                 end
%             end
%             figure;imshow(Cl(:,:,2),[0 1]);hold on;
%             for i=1:length(path_vol) 
% 
%                 if eq(path_vol(i),2)
%                     idx2=idx(i);
%                     plot3(X(idx2),Y(idx2),Z(idx2),'-','Linewidth',2, 'Color',xroma);
%                 end
%             end
%             figure;imshow(Cl(:,:,3),[0 1]);hold on;
%             for i=1:length(path_vol) 
% 
%                 if eq(path_vol(i),3)
%                     idx3=idx(i);
%                     plot3(X(idx3),Y(idx3),Z(idx3),'-','Linewidth',2, 'Color',xroma);
%                 end
%             end
% 
%             figure;imshow(Cl(:,:,4),[0 1]);hold on;
%             for i=1:length(path_vol) 
% 
%                 if eq(path_vol(i),4)
%                     idx4=idx(i);
%                     plot3(X(idx4),Y(idx4),Z(idx4),'-','Linewidth',2, 'Color',xroma);
%                 end
%             end
% 
%     end
    
end