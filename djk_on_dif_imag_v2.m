%efarmogi tou djkstra (2D) se magnitiki me diffussion 
%Epilogi ton simeion me xrisi tou getpts

Slice=15;
A0=vol_3D(:,:,15); 

index_cl=find(A0>150 & cl(:,:,Slice)>0.3);
[x0,y0]=ind2sub(size(cl),index_cl);
fig1=figure;
imshow(cl(:,:,Slice),[0,1]); title('stis 2 diastaseis');
hold on; 

A1=Diffusion_tensor_eigenvectors(:,:,Slice,1);
A2=Diffusion_tensor_eigenvectors(:,:,Slice,2);

quiver(y0,x0,A1(index_cl),A2(index_cl));

[x,y]=getpts(fig1); 

% x=100;
% y=160;
% x_t=110;
% y_t=120;
x1=round(x(1));
y1=round(y(1));
x2=round(x(end));
y2=round(y(end));
% x1=round(x);
% y1=round(y);
% x2=round(x_t);
% y2=round(y_t);


%metatropi ton dianismaton apo cartesianes se polikes :

%angle_th_eigvctr einai i gonia se rad
%radius_eigvctr to mikos t?? dianismatoa
[angle_th_eigvctr,radius_eigvctr]=cart2pol(Diffusion_tensor_eigenvectors(:,:,Slice,1),...
                                            Diffusion_tensor_eigenvectors(:,:,Slice,2));

[d,px,py,QS]=djk_fun_on_dif_imag1(A0,y1,x1,Diffusion_tensor_eigenvectors(:,:,Slice,:),...
                                        angle_th_eigvctr,radius_eigvctr,cl(:,:,Slice));
[pathx,pathy,b]=backtrack_fun_vers2(d,A0,px,py,y2,x2);


imshow(cl(:,:,Slice),[0,1]); title('path');
hold on;
quiver(y0,x0,A1(index_cl),A2(index_cl));
plot(pathy,pathx,'r-','Linewidth',2);
plot(pathy,pathx,'r-');

