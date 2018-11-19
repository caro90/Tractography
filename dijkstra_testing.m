%Select a slice
slice1=20;

f1=fib.dimension(1);
f2=fib.dimension(2);
[X,Y]=meshgrid(1:f1,1:f2);
Vx=fib.dir0(:,:,slice1,1);
Vy=fib.dir0(:,:,slice1,2);
fa=fib.fa0(:,:,slice1);

%A0=vol_3D(:,:,slice1); 
%%Threshold for quiver
%index_fa=find(A0>150 & fa>0.1);
%[x0,y0]=ind2sub(size(fa),index_fa);
%fig1=figure1;
%figure;imshow(fib.fa0(:,:,20))
%hold on; 
%A1=vars.dir0(:,:,slice1,1);
% A2=vars.dir0(:,:,slice1,2);
%quiver(y0,x0,A1(index_cl),A2(index_cl));
%[x_path,y_path]=getpts(figure(1)); 

for i=1:3
    %Choosing a point from the dataset
    x1=trk.tracts(1,i);
    y1=trk.tracts(2,i);
    x2=trk.tracts( 1,trk.length(i));
    y2=trk.tracts( 2,trk.length(i));
    
    x_path=[x1,x2]; 
    y_path=[y1,y2]; 
    i1=round(x_path(1));
    j1=round(y_path(1));
    i2=round(x_path(end));
    j2=round(y_path(end));
    % Dijkstra 
    [d(:,:,i),pi(:,:,i),pj(:,:,i),QS1]=dijkstra_function(j1,i1,X,Y,Vx,Vy,fa);

    % plotting and backtracking 
    h1=plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',10);
    set(h1,'MarkerFaceColor',[0,1,0]);
    cmap=colormap(colorcube);

    % imshow(fa,[0,1]); title('path');
    % hold on;
    % quiver(y0,x0,A1(index_cl),A2(index_cl));
    % plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',15);

    % Plotting the path
     for k=1:1   
        xroma=cmap(mod(k+6,size(cmap,1))+1,:);
        [path_lin,path_col,B_backtrack]=dijkstra_backtrack_fun(d,pi,pj,j2,i2); 
        idx=sub2ind(size(X),path_lin,path_col); 
        hold on;
        plot(X(j2,i2),Y(j2,i2),'p','MarkerSize',15,'Color',xroma);
        h=plot(X(idx),Y(idx),'-','Linewidth',3, 'Color',xroma); 
        drawnow;
     end
     
%     path_lin=fliplr(path_lin');
%     path_col=fliplr(path_col');
%     %Error calculation
%     for j=1: length(path_lin)
%         error(j,1)=path_lin(j)-idxRounded(j,1);
%         error(j,2)=path_col(j)-idxRounded(j,2);
%     end
    clear path_lin path_col
end
clear slice1 f1 f2 X Y Anisotropy Vx Vy fa x_path y_path i1 j1 i2 j2 angle_TH Radius k cmap xroma x1 y1 x2 y2 i
    
