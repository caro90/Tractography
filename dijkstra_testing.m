%%% Efarmogi tou dijkstra se VectorFields 2D kai real data
% vol_3D=img_double;
% f1=size(Diffusion_tensor,1);
% f2=size(Diffusion_tensor,2);

%% ----Real data--------
% slice1=37;
% %cl normalized:
% Cl=cl_normal(:,:,slice1);
% 
% Vx=Diffusion_tensor_eigenvectors(:,:,slice1,1);
% Vy=Diffusion_tensor_eigenvectors(:,:,slice1,2);
% [X,Y]=meshgrid(1:f1,1:f2);%f2:-1:1);
% A0=vol_3D(:,:,slice1); 
% % 
% %% threshold gia na min emfanisei ola ta idiodianismata sto quiver
% index_cl=find(A0>150 & Cl>0.1);
%[x0,y0]=ind2sub(size(Cl),index_cl);
%fig1=figure1;
%imshow(Cl,[0 1]);
hold on; 
% A1=Diffusion_tensor_eigenvectors(:,:,slice1,1);
% A2=Diffusion_tensor_eigenvectors(:,:,slice1,2);
% quiver(y0,x0,A1(index_cl),A2(index_cl));

[x_path,y_path]=getpts(figure(1) ); 
%----------------------

% Epilogi simeiou gia Dijkstra xoris getpts
%  x_path=[59,5]; 
%  y_path=[5,60]; 

%%
i1=round(x_path(1));
j1=round(y_path(1));
i2=round(x_path(end));
j2=round(y_path(end));

Vx=V1x;
Vy=V1y;
Vx=reshape(Vx,N,N);
Vy=reshape(Vy,N,N);

%%
L1max=max(L1(:));

%%
[angle_TH,Radius] = cart2pol(Vx,Vy);

%% dijkstra function call
[d,pi,pj,QS1]=dijkstra_function(j1,i1,X,Y,Vx,Vy,Cl,Cp,L1max,L1);%,trace(:,:,slice1),FA(:,:,slice1),MD(:,:,slice1),vol_3D(:,:,slice1),Cp); 

%% plotting kai backtrack function call
h1=plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',10);
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);

% i dests kai i for xrisimopoieitai gia na kanoume polles fores
% backtrack apo diafora simeia meta ton dijkstra.
% epilegoume polla tyxaia simeia stin eikona:
% dests=randint(50,2,[1+2,size(X,1)-2]);
% dests=[130 20 ; 140 40 ; 20 100 ; 22 150 ; 19 155 ; ];%; 11 101 ; 10 90 ; 12 98  ; 13 95 ; 15 80];
%dests=[ 100 100;  90 101 ;  95 98 ;  100  90 ;  99 98 ;  85 101; 80 90];

% Gia tin efarmogi se real data (oxi vectorfields)
%     imshow(Cl,[0,1]); title('path');
%     hold on;
%     quiver(y0,x0,A1(index_cl),A2(index_cl));
%     plot(X(j1,i1),Y(j1,i1),'o','MarkerSize',15);

 for k=1:1%length(dests(:,2))   
    xroma=cmap(mod(k+6,size(cmap,1))+1,:);
%     i2=dests(k,1);
%     j2=dests(k,2);
      [path_lin,path_col,B_backtrack]=dijkstra_backtrack_fun(d,pi,pj,j2,i2); 
  
    
    %% plotting to path
    idx=sub2ind(size(X),path_lin,path_col); 
    hold on;
    plot(X(j2,i2),Y(j2,i2),'p','MarkerSize',15,'Color',xroma);
    h=plot(X(idx),Y(idx),'-','Linewidth',3, 'Color',xroma); 
    drawnow;
    
 end
    
    
