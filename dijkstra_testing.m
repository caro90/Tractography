slice1=20;

f1=fib.dimension(1);
f2=fib.dimension(2);
[X,Y]=meshgrid(1:f1,1:f2);
Anisotropy=fib.fa0(:,:,slice1);
Vx=fib.dir0(:,:,slice1,1);
Vy=fib.dir0(:,:,slice1,2);
fa=fib.fa0(:,:,slice1);

%A0=vol_3D(:,:,slice1); 
%threshold gia na min emfanisei ola ta idiodianismata sto quiver
%index_cl=find(A0>150 & Cl>0.1);
%[x0,y0]=ind2sub(size(Cl),index_cl);
%fig1=figure1;
%imshow(Cl,[0 1]);
%hold on; 
% A1=vars.dir0(:,:,slice1,1);
% A2=vars.dir0(:,:,slice1,2);
% quiver(y0,x0,A1(index_cl),A2(index_cl));
%[x_path,y_path]=getpts(figure(1) ); 

% Manually choose a point
x_path=[49,5]; 
y_path=[54,75]; 
i1=round(x_path(1));
j1=round(y_path(1));
i2=round(x_path(end));
j2=round(y_path(end));

%??
[angle_TH,Radius] = cart2pol(Vx,Vy);

% dijkstra 
[d,pi,pj,QS1]=dijkstra_function(j1,i1,X,Y,Vx,Vy,fa);

% plotting and backtracking 
h1=plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',10);
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);

%Use dest to backtrack from many points
% dests=randint(50,2,[1+2,size(X,1)-2]);
% dests=[130 20 ; 140 40 ; 20 100 ; 22 150 ; 19 155 ; ];

% imshow(fa,[0,1]); title('path');
% hold on;
% quiver(y0,x0,A1(index_cl),A2(index_cl));
% plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',15);

 for k=1:1%length(dests(:,2))   
    xroma=cmap(mod(k+6,size(cmap,1))+1,:);
%     i2=dests(k,1); %j2=dests(k,2);
    [path_lin,path_col,B_backtrack]=dijkstra_backtrack_fun(d,pi,pj,j2,i2); 
    % Plotting the path
    idx=sub2ind(size(X),path_lin,path_col); 
    hold on;
    plot(X(j2,i2),Y(j2,i2),'p','MarkerSize',15,'Color',xroma);
    h=plot(X(idx),Y(idx),'-','Linewidth',3, 'Color',xroma); 
    drawnow;
 end
    
    
