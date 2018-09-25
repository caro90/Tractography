%%% Efarmogi tou dijkstra se VectorFields 2D kai real data

% 
% if  -1>0
%     %% Dimiourgia ton Vector fields
%     [X,Y]=meshgrid(-2:0.02:2,-2:0.02:2);
%     Vx=-Y;
%     Vy=X;
% 
%     Vx=sin(X);
%     Vy=sin(2*Y);
%     [Vx,Vy]=gradient(peaks(X,Y));   %HELP
%     Cl=sqrt(Vx.^2+Vy.^2);
%     Cl=Cl/max(Cl(:));
% else
%     %% Real data
%     
%     %cl normalized:
%     temp_cl_new=(Diffusion_tensor_eigenvalues_sorted(:,:,16,1)-Diffusion_tensor_eigenvalues_sorted(:,:,16,2))./ ...
%         sqrt(Diffusion_tensor_eigenvalues_sorted(:,:,16,1).^2+Diffusion_tensor_eigenvalues_sorted(:,:,16,2).^2+...
%         Diffusion_tensor_eigenvalues_sorted(:,:,16,3).^2);
%     temp_cl_new(isnan(temp_cl_new))=0;
%     Vx=Diffusion_tensor_eigenvectors(:,:,16,1);
%     Vy=Diffusion_tensor_eigenvectors(:,:,16,2);
%     [X,Y]=meshgrid(1:256,1:256);
%     Cl=temp_cl_new;
%     
% 
%         
% end

%fig1=figure;
% imshow(cl(:,:,16),[ ]);
hold on;quiver(X,Y,Vx,Vy);
%[x_path,y_path]=getpts(fig1); 

 x_path=[35,4]; 
 y_path=[38,37];
 

i1=round(x_path(1));
j1=round(y_path(1));
i2=round(x_path(end));
j2=round(y_path(end));

%% metatropi ton dianismaton apo cartesianes se polikes :
%angle_th_eigvctr einai i gonia se rad
%radius_eigvctr to mikos ton dianismaton
[angle_th_eigvctr,radius_eigvctr]=cart2pol(Vx,Vy);
axis equal

%% dijkstra calling
[d,pi,pj,QS]=dijkstra_testing_with_VectorsFields_djk_fun(i1,j1,X,Y,Vx,Vy,Cl); 

%% plotting kai backtrack calling
h1=plot(X(i1,j1),Y(i1,j1),'o','MarkerSize',10);
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);
%tyxaia simeia:
dests=randint(200,2,[1+2,size(X,1)-2]);

for k=1:1     
    xroma=cmap(mod(k+1,size(cmap,1))+1,:);
    %i2=dests(k,1);
    %j2=dests(k,2);
                                                   
    [path_lin,path_col,B_backtrack]=dijkstra_testing_with_VectorsFields_backtrack_fun(d,pi,pj,i2,j2); 
    idx=sub2ind(size(X),path_lin,path_col);    
    plot(X(i2,j2),Y(i2,j2),'p','MarkerSize',5,'Color',xroma);
    h=plot(X(idx),Y(idx),'-','Linewidth',2, 'Color',xroma); 
    drawnow;
end