%2D dijkstra function

%Diorthosi tou kodika djk_fun_on_dif_imag.m :

function [d,px,py,QS]=djk_fun_on_dif_imag1(I,x0,y0,D_temp_eigenvectors,angle_th_eigvctr,radius_eigvctr,cl)
D_temp_eigenvectors=reshape(D_temp_eigenvectors,256,256,9);

a=double(I);
[N,M]=size(a);
v0=I(x0,y0);
k=0;

d=ones(size(a))*1e9;
px=-ones(size(a));
py=-ones(size(a));
d(x0,y0)=0;

px(x0,y0)=x0;
py(x0,y0)=y0;

QS=ones(size(a));
m0=max(max(d));

tic;
while(sum(sum(QS))>0)
    
    
%     m=m0;       
%     for i=1:256
%         for j=1:256
%             if QS(i,j)==1 && d(i,j)<=m 
%                 xmin=i;
%                 ymin=j;
%                 m=d(i,j);
%             end
%         end
%     end
%     QS(xmin,ymin)=0;
    
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));
    QS(idx1(ii))=0;
    [xmin,ymin]=ind2sub([256,256],idx1(ii));

%-----------------------------------------------------------------------------------------------  
f=angle_th_eigvctr(xmin,ymin);

    a_vector= [D_temp_eigenvectors(xmin,ymin,1) D_temp_eigenvectors(xmin,ymin,2) ] ;
    metro_a=radius_eigvctr(xmin,ymin);        

  
    %Katanomi ton varon metavasis sta 8 geitonika pixel :
 
    %1
    if ( xmin-1 > 0) &&  ymin+1<=M  
        
        %to b_vector einai to idiodianisma sto epomeno pixel
        b_vector=[ D_temp_eigenvectors(xmin-1,ymin+1,1) D_temp_eigenvectors(xmin-1,ymin+1,2) ];
        esot_gin=dot(a_vector,b_vector); 
        metro_b=radius_eigvctr(xmin-1,ymin+1); 
        cos_th=(esot_gin/( metro_a * metro_b ));        
        %v_unit=[cos(TH),sin(Th)];
        % cost2=acos(sum(a_vector.*v_unit));
        
        TH=-pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;% gia na doume poso omoies einai ooi th kai f
        dd=(1-cos_th)*cost2;
        if d(xmin-1,ymin+1)>d(xmin,ymin) +dd
%             d (xmin-1,ymin+1)=d(xmin,ymin)+(1-cos_th);%*cost2;%*(1-cl(xmin-1,ymin+1));%*cost2;
            d (xmin-1,ymin+1)=d(xmin,ymin)+dd;
            px(xmin-1,ymin+1)=xmin;
            py(xmin-1,ymin+1)=ymin;
        end
    end
    
    %2
    if  ymin+1<=M         
        b_vector=[ D_temp_eigenvectors(xmin,ymin+1,1) D_temp_eigenvectors(xmin,ymin+1,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin+1); 
        cos_th=(esot_gin/( metro_a * metro_b ));                
        TH=0;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5; %allagi dn eixe *1.5
        dd=(1-cos_th)*cost2; %*(1-cl(xmin,ymin+1));%*cost2;
        if d(xmin,ymin+1)>d(xmin,ymin) + dd
            d(xmin,ymin+1)=d(xmin,ymin) + dd;
            px(xmin,ymin+1)=xmin;
            py(xmin,ymin+1)=ymin;
        end
    end
    %3
    
     if  xmin+1<=N && ymin+1<=M          
        b_vector=[ D_temp_eigenvectors(xmin+1,ymin+1,1) D_temp_eigenvectors(xmin+1,ymin+1,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin+1,ymin+1);    
        cos_th=(esot_gin/( metro_a * metro_b ));
 
        TH=pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        dd=(1-cos_th)*cost2;
        if d(xmin+1,ymin+1)>d(xmin,ymin)+dd
            d(xmin+1,ymin+1)=d(xmin,ymin)+dd;
            px(xmin+1,ymin+1)=xmin;
            py(xmin+1,ymin+1)=ymin;
        end
     end
     
     %4
     if  xmin+1<=N &&  ymin-1>0          
        b_vector=[D_temp_eigenvectors(xmin+1,ymin-1,1),D_temp_eigenvectors(xmin+1,ymin-1,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin+1,ymin-1);   
        cos_th=(esot_gin/( metro_a * metro_b ));
    
        TH=3*pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        dd=(1-cos_th)*cost2;
        if d(xmin+1,ymin-1)>d(xmin,ymin) +dd
            d(xmin+1,ymin-1)=d(xmin,ymin)+dd;  %  (1-cos_th);%*cost2;%*(1-cl(xmin+1,ymin-1));%*cost2;
            px(xmin+1,ymin-1)=xmin;
            py(xmin+1,ymin-1)=ymin;
        end
     end
     
     %5
     if xmin+1<=N        
         b_vector=[D_temp_eigenvectors(xmin+1,ymin,1),D_temp_eigenvectors(xmin+1,ymin,2) ];
         esot_gin=dot(a_vector,b_vector);
         metro_b=radius_eigvctr(xmin+1,ymin);
         cos_th=(esot_gin/( metro_a * metro_b ));
         
         TH=pi/2;
         cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
         dd=(1-cos_th)*cost2;
         if d(xmin+1,ymin)>d(xmin,ymin) +dd
             d(xmin+1,ymin)=d(xmin,ymin)+dd;  %  (1-cos_th);%*cost2;%*(1-cl(xmin+1,ymin));%*cost2;
             px(xmin+1,ymin)=xmin;
             py(xmin+1,ymin)=ymin;
         end
     end
     
     %6
    if xmin-1>0 &&  ymin-1>0   
        b_vector=[D_temp_eigenvectors(xmin-1,ymin-1,1),D_temp_eigenvectors(xmin-1,ymin-1,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin-1,ymin-1);  
        cos_th=(esot_gin/( metro_a * metro_b ));
    
        TH=-3*pi/4;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5;
        dd=(1-cos_th)*cost2;
        if d(xmin-1,ymin-1)>d(xmin,ymin)+dd
            d (xmin-1,ymin-1)=d(xmin,ymin)+dd;  %  (1-cos_th);%*cost2;%*(1-cl(xmin-1,ymin-1));%*cost2;
            px(xmin-1,ymin-1)=xmin;
            py(xmin-1,ymin-1)=ymin;
        end
    end
    
    %7
    if ymin-1>0         
        b_vector=[D_temp_eigenvectors(xmin,ymin-1,1),D_temp_eigenvectors(xmin,ymin-1,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin,ymin-1);    
        cos_th=(esot_gin/( metro_a * metro_b ));
        
        TH=-pi;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5; 
        dd=(1-cos_th)*cost2;
        if d(xmin,ymin-1)>d(xmin,ymin)+dd
            d(xmin,ymin-1)=d(xmin,ymin)+dd;  %  (1-cos_th);%*cost2;%*(1-cl(xmin,ymin-1));%*cost2;
            px(xmin,ymin-1)=xmin;
            py(xmin,ymin-1)=ymin;
        end
    end
    
    %8
    if xmin-1>0      
        b_vector=[D_temp_eigenvectors(xmin-1,ymin,1),D_temp_eigenvectors(xmin-1,ymin,2) ];
        esot_gin=dot(a_vector,b_vector);
        metro_b=radius_eigvctr(xmin-1,ymin);    
        cos_th=(esot_gin/( metro_a * metro_b ));
        
        TH=-pi/2;
        cost2=[abs((cos(f)-cos(TH))).*abs((sin(f)-sin(TH)))]*1.5; 
        dd=(1-cos_th)*cost2;
        if d(xmin-1,ymin)>d(xmin,ymin) +dd
            d(xmin-1,ymin)=d(xmin,ymin)+(1-cos_th);%*cost2;%*(1-cl(xmin-1,ymin));%*cost2;
            px(xmin-1,ymin)=xmin;
            py(xmin-1,ymin)=ymin;
        end
    end
   
    k=k+1;
    if mod(k,100)==0
        k;
    end
end
  toc  
end