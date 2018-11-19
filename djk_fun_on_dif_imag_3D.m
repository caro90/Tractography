%3D dijkstra function

function [d,px,py,pz,QS]=djk_fun_on_dif_imag_3D(I,x0,y0,z0,Diffusion_tensor_eigenvectors,...
                                                           angle_th_eigvctr,radius_eigvctr,cl)

%to I einai 3D
a=double(I);
v0=I(x0,y0,z0);

d=ones(size(a))*1e9;
px=-ones(size(a));
py=-ones(size(a));
%epilogi arxikou simeiou:
d(x0,y0,z0)=0; 

px(x0,y0,z0)=x0;
py(x0,y0,z0)=y0;
pz(x0,y0,z0)=z0;

QS=ones(size(a));
m0=max(max(max(d)));

tic;
while sum((sum(sum(QS)))) >0
    m=m0;
    for i=180:256
        for j=180:256
            for k=1:2  %gia ton z aksona // prosorina gia 2 tomes
                if QS(i,j,k)==1 && d(i,j,k)<=m 
                    xmin=i;
                    ymin=j;
                    zmin=k;
                    
                    m=d(i,j,k);
                end
            end
        end
    end
 QS(xmin,ymin,zmin)=0; %mark processed pixel
 
 
%-----------------------------------------------------------------------------------------------  
f=angle_th_eigvctr(xmin,ymin,zmin);
a_vector= [Diffusion_tensor_eigenvectors(xmin,ymin,zmin,1) Diffusion_tensor_eigenvectors(xmin,ymin,zmin,2) ...
                                                           Diffusion_tensor_eigenvectors(xmin,ymin,zmin,3) ] ;
metro_a=radius_eigvctr(xmin,ymin,zmin);    
 %1
    if (( xmin-1>0) &&  ymin+1<=size(a,2) && d(xmin-1,ymin+1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin-1,ymin+1,zmin)) )
        
        b_vector=[ Diffusion_tensor_eigenvectors(xmin-1,ymin+1,1) Diffusion_tensor_eigenvectors(xmin-1,ymin+1,2)...
                                                                  Diffusion_tensor_eigenvectors(xmin-1,ymin+1,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin-1,ymin+1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));
        
        %cos_th=round(cos_th);
      
        TH=-pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        d (xmin-1,ymin+1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);    %*cost2*(1-cl(xmin-1,ymin+1,zmin));
        px(xmin-1,ymin+1,zmin)=xmin;
        py(xmin-1,ymin+1,zmin)=ymin;
        pz(xmin-1,ymin+1,zmin)=zmin;

    end
 %2
     if ymin+1<=size(a,2) && d(xmin,ymin+1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin,ymin+1,zmin)) 

        b_vector=[ Diffusion_tensor_eigenvectors(xmin,ymin+1,zmin,1) Diffusion_tensor_eigenvectors(xmin,ymin+1,zmin,2)...
                                                                    Diffusion_tensor_eigenvectors(xmin,ymin+1,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin+1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));
        
        %cos_th=round(cos_th);
        
        TH=0;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d(xmin,ymin+1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);   %*cost2*(1-cl(xmin,ymin+1,zmin));
        px(xmin,ymin+1,zmin)=xmin;
        py(xmin,ymin+1,zmin)=ymin;
        pz(xmin,ymin+1,zmin)=zmin;
    end
   %3 
     if ( xmin+1<=size(a,1) &&  ymin+1<=size(a,2) && d(xmin+1,ymin+1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin+1,ymin+1,zmin)) )
        
        b_vector=[Diffusion_tensor_eigenvectors(xmin+1,ymin+1,zmin,1) Diffusion_tensor_eigenvectors(xmin+1,ymin+1,zmin,2)...
                                                                       Diffusion_tensor_eigenvectors(xmin+1,ymin+1,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin+1,ymin+1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
    
        TH=pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        d (xmin+1,ymin+1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);    %*cost2*(1-cl(xmin+1,ymin+1,zmin);
        px(xmin+1,ymin+1,zmin)=xmin;
        py(xmin+1,ymin+1,zmin)=ymin;
        pz(xmin+1,ymin+1,zmin)=zmin;

     end
    %4
    if  xmin+1<=size(a,1) &&  ymin-1>0 && d(xmin+1,ymin-1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin+1,ymin-1,zmin))
      
        b_vector=[Diffusion_tensor_eigenvectors(xmin+1,ymin-1,zmin,1) Diffusion_tensor_eigenvectors(xmin+1,ymin-1,zmin,2)...
                                                                Diffusion_tensor_eigenvectors(xmin+1,ymin-1,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin+1,ymin-1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
    
        TH=3*pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        d (xmin+1,ymin-1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);    %*cost2*(1-cl(xmin+1,ymin-1,zmin));
        px(xmin+1,ymin-1,zmin)=xmin;
        py(xmin+1,ymin-1,zmin)=ymin;
        pz(xmin+1,ymin-1,zmin)=zmin;

    end
    %5
    if xmin+1<=size(a,1) && d(xmin+1,ymin,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin+1,ymin,zmin))
        
        b_vector=[ Diffusion_tensor_eigenvectors(xmin+1,ymin,zmin,1) Diffusion_tensor_eigenvectors(xmin+1,ymin,zmin,2)...
                                                                    Diffusion_tensor_eigenvectors(xmin+1,ymin,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin+1,ymin,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
        
        TH=pi/2;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d (xmin+1,ymin,zmin)=d(xmin,ymin,zmin)+(1-cos_th);      %*cost2*(1-cl(xmin+1,ymin,zmin));    
        px(xmin+1,ymin,zmin)=xmin;
        py(xmin+1,ymin,zmin)=ymin;
        pz(xmin+1,ymin,zmin)=zmin; 
    end
    %6
    if( xmin-1>0 &&  ymin-1>0 && d(xmin-1,ymin-1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin-1,ymin-1,zmin)) )
        
        b_vector=[ Diffusion_tensor_eigenvectors(xmin-1,ymin-1,zmin,1) Diffusion_tensor_eigenvectors(xmin-1,ymin-1,zmin,2)...
                                                                Diffusion_tensor_eigenvectors(xmin-1,ymin-1,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin-1,ymin-1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
    
        TH=-3*pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        d(xmin-1,ymin-1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);     %*cost2*(1-cl(xmin-1,ymin-1,zmin));
        px(xmin-1,ymin-1,zmin)=xmin;
        py(xmin-1,ymin-1,zmin)=ymin;
        pz(xmin-1,ymin-1,zmin)=zmin;
    end
    %7
    if ymin-1>0 && d(xmin,ymin-1,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin,ymin-1,zmin))
        
        b_vector=[ Diffusion_tensor_eigenvectors(xmin,ymin-1,zmin,1) Diffusion_tensor_eigenvectors(xmin,ymin-1,zmin,2)... 
                                                                 Diffusion_tensor_eigenvectors(xmin,ymin-1,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin-1,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));
 
        %cos_th=round(cos_th);
        
        TH=-pi;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d(xmin,ymin-1,zmin)=d(xmin,ymin,zmin)+(1-cos_th);       %*cost2*(1-cl(xmin,ymin-1,zmin));
        px(xmin,ymin-1,zmin)=xmin;
        py(xmin,ymin-1,zmin)=ymin;
        pz(xmin,ymin-1,zmin)=zmin;
    end
    %8
    if xmin-1>0 && d(xmin-1,ymin,zmin)>d(xmin,ymin,zmin)+abs(v0-a(xmin-1,ymin,zmin))
      
        b_vector=[ Diffusion_tensor_eigenvectors(xmin-1,ymin,zmin,1) Diffusion_tensor_eigenvectors(xmin-1,ymin,zmin,2)...
                                                                  Diffusion_tensor_eigenvectors(xmin-1,ymin,zmin,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin-1,ymin,zmin);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
        
        TH=-pi/2;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d(xmin-1,ymin,zmin)=d(xmin,ymin,zmin)+(1-cos_th);       %*cost2*(1-cl(xmin-1,ymin,zmin));
        px(xmin-1,ymin,zmin)=xmin;
        py(xmin-1,ymin,zmin)=ymin;
        pz(xmin-1,ymin,zmin)=zmin;
    end
    
    %--------gia na kinithei pano kato ston z aksona-------------
    %allagi zmin+1<=27
    %9
    if zmin+1<=2 && d(xmin,ymin,zmin+1)>d(xmin,ymin,zmin)+abs(v0-a(xmin,ymin,zmin+1))  
       
        b_vector=[ Diffusion_tensor_eigenvectors(xmin,ymin,zmin+1,1) Diffusion_tensor_eigenvectors(xmin,ymin,zmin+1,2)...
                                                                Diffusion_tensor_eigenvectors(xmin,ymin,zmin+1,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin,zmin+1);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
        
        TH=-pi/2;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d(xmin,ymin,zmin+1)=d(xmin,ymin,zmin)+(1-cos_th);       %*cost2*(1-cl(xmin,ymin,zmin+1));
        px(xmin,ymin,zmin+1)=xmin;
        py(xmin,ymin,zmin+1)=ymin;
        pz(xmin,ymin,zmin+1)=zmin;
    end
    %10
    if zmin-1>0 && d(xmin,ymin,zmin-1)>d(xmin,ymin,zmin)+abs(v0-a(xmin,ymin,zmin-1))  
       
        b_vector=[ Diffusion_tensor_eigenvectors(xmin,ymin,zmin-1,1) Diffusion_tensor_eigenvectors(xmin,ymin,zmin-1,2)...
                                                                Diffusion_tensor_eigenvectors(xmin,ymin,zmin-1,3) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin,zmin-1);
        cos_th=abs(esot_gin/( metro_a * metro_b ));

        %cos_th=round(cos_th);
        
        TH=-pi/2;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))];
        d(xmin,ymin,zmin-1)=d(xmin,ymin,zmin)+(1-cos_th);       %*cost2*(1-cl(xmin,ymin,zmin-1));
        px(xmin,ymin,zmin-1)=xmin;
        py(xmin,ymin,zmin-1)=ymin;
        pz(xmin,ymin,zmin-1)=zmin;
    end
    
    %--------gia na kinithei pano kato ston z-------------
 
end
toc;
sum((sum(sum(QS))))

end