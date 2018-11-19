%Efarmogi tou dijkstra 3D

% Slice=7;
% A0=vol_3D(:,:,40); 
% 
% index_cl=find(A0>150 & cl(:,:,Slice)>0.4);
% [x0,y0]=ind2sub(size(cl),index_cl);
% fig1=figure;
% imshow(cl(:,:,Slice),[0,1]); title('stis 3 diastaseis');
% hold on; 
% 
% A1=Diffusion_tensor_eigenvectors(:,:,Slice,1);
% A2=Diffusion_tensor_eigenvectors(:,:,Slice,2);
% A3=Diffusion_tensor_eigenvectors(:,:,Slice,3);
% quiver(y0,x0,A1(index_cl),A2(index_cl),);



%angle_th_eigvctr : oi gonies gia ola ta pixel olon ton tomon se polikes sintetagmenes
%radius_eigvctr   : to metro tous.
[angle_th_eigvctr,radius_eigvctr]=cart2pol(Diffusion_tensor_eigenvectors(:,:,:,1), ...
                                    Diffusion_tensor_eigenvectors(:,:,:,2),Diffusion_tensor_eigenvectors(:,:,:,3)); 

%x0 ,y0 , z0 to tyxaio arxiko simeio apo opou thelo na ksekinisei o dijkstra
x0=190;y0=200;z0=1;
%x2, y2 , z2 einai to tyxaio teliko simeio
x2=250;y2=240;z2=2;

%To S0 einai oles oi tomes gia gradient=0. 
%prosorini allagi (den exoume dosei olo to dataset alla mono 2 tomes):
new=S0(:,:,1:2);

[d,px,py,pz,QS]=djk_fun_on_dif_imag_3D(new,y0,x0,z0,Diffusion_tensor_eigenvectors,angle_th_eigvctr,radius_eigvctr,cl);

[pathx,pathy,pathz,b]=backtrack_fun_3D(new,px,py,pz,y2,x2,z2);

