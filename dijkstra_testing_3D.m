%% Efarmogi tou dijkstra se VectorFields 3D kai real data
%Version 11/2/17
% flag==0 for synthetic data
% flag==1 for real data

flag=1;

% gia sigekrimenou typou dataset 
% vol_3D=img_double;
% vol_3D=img_double(:,:,:,1);

if (flag==1 )
%% fortosi idiodianismaton gia real data
% slice1 kai slice 2 einai oi tomes metaksi ton 
% opoion theloume na treksei o dijkstra
        slice1=44;
        slice2=50;
        Vx=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,1);
        Vy=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,2);
        Vz=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,3);

        Vx3=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,1);
        Vy3=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,2);
        Vz3=Diffusion_tensor_eigenvectors(:,:,slice1:slice2,3);
        
        f1=size(Diffusion_tensor,1);
        f2=size(Diffusion_tensor,2);
        f3=size(Diffusion_tensor,3);
           
        % f1 f2 reverse
        [X,Y,Z]=meshgrid(1:f2,1:f1,slice1:slice2);
        Cl=cl_normal(:,:,slice1:slice2);
        Cp=cp_normal(:,:,slice1:slice2);
        Cs=cs_normal(:,:,slice1:slice2);
        Cs(Cs<0)=0;
        
end

%% epilogi simeiou  
x_path=[40,38]; 
y_path=[65,34];
% to z setarete os sinolikos arithmos ton tomon
z_path=[1,7];

i1=round(x_path(1));
j1=round(y_path(1));
z1=round(z_path(1));
i2=round(x_path(end));
j2=round(y_path(end));
z2=round(z_path(end));


SliceThickness=1;
[d,p_i,pj,pz,QS]=dijkstra_function_3D_version2(i1,j1,z1,X,Y,Z,Vx,Vy,Vz,Vx3,Vy3,Vz3,...
                            Cl,Cp,Cs,FA(:,:,slice1:slice2),trace(:,:,slice1:slice2),...
                   MD(:,:,slice1:slice2),vol_3D(:,:,slice1:slice2),SliceThickness,flag);


%% plotting kai backtrack function call
h1=plot3(X(i1,j1,z1),Y(i1,j1,z1),Z(i1,j1,z1),'o','MarkerSize',10);
set(h1,'MarkerFaceColor',[0,1,0]);
cmap=colormap(colorcube);
% i dests kai i for xrisimopoieitai gia na kanoume polles fores
% backtrack apo diafora simeia meta ton dijkstra.
% epilegoume polla tyxaia simeia stin eikona:
dests=randint(100,2,[1+2,size(X,1)-2]);

for k=1:1      
    xroma=cmap(mod(k+1,size(cmap,1))+1,:);
    % i2=dests(k,1);
    % j2=dests(k,2);
                                                   
    [path_lin,path_col,path_vol,B_backtrack]=dijkstra_backtrack_fun_3D(d,p_i,pj,pz,i2,j2,z2); 
    idx=sub2ind(size(X),path_lin,path_col,path_vol);    
    plot3(X(i2,j2,z2),Y(i2,j2,z2),Z(i2,j2,z2),'p','MarkerSize',10,'Color',xroma);
    %ayto to plot xrisimopoieitai gia synthetic 3D data:
    %if (flag==0)
        figure;
        h=plot3(X(idx),Y(idx),Z(idx),'-','Linewidth',2, 'Color',xroma);
        drawnow;
    %end
    
    %% ---------------- plotting test gia real data -----------
    if (flag==3)
            figure;imshow(Cl(:,:,1),[0 1]);hold on;
            for i=1:length(path_vol) 

                if eq(path_vol(i),1)
                    idx1=idx(i);
                    plot3(X(idx1),Y(idx1),Z(idx1),'-','Linewidth',2, 'Color',xroma);
                end
            end
            figure;imshow(Cl(:,:,2),[0 1]);hold on;
            for i=1:length(path_vol) 

                if eq(path_vol(i),2)
                    idx2=idx(i);
                    plot3(X(idx2),Y(idx2),Z(idx2),'-','Linewidth',2, 'Color',xroma);
                end
            end
            figure;imshow(Cl(:,:,3),[0 1]);hold on;
            for i=1:length(path_vol) 

                if eq(path_vol(i),3)
                    idx3=idx(i);
                    plot3(X(idx3),Y(idx3),Z(idx3),'-','Linewidth',2, 'Color',xroma);
                end
            end

            figure;imshow(Cl(:,:,4),[0 1]);hold on;
            for i=1:length(path_vol) 

                if eq(path_vol(i),4)
                    idx4=idx(i);
                    plot3(X(idx4),Y(idx4),Z(idx4),'-','Linewidth',2, 'Color',xroma);
                end
            end

    end
    
end