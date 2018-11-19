%2D Dijkstra
function [d,pi,pj,QS1]=dijkstra_function(i0,j0,X,Y,Vx,Vy,fa)
global p0;
global f1;
global f2;
global f3;

%cost table:
d=ones(size(Vx))*1e9;
d(i0,j0)=0;
%pi kai pj (previous) keeps the previous pixel(position index) 
pi=-ones(size(Vx));
pj=-ones(size(Vx));
%Initialize the first point
pi(i0,j0)=i0;
pj(i0,j0)=j0;

%In Qs we set 0 to the pixel we processed 
QS=ones(size(Vx));
[N,M]=size(pi);
% Trace,FA and MD thresholds to restrict dijkstra inside the region of the head
% for i=1:N
%     for j=1:M
%         if vol_3D(i,j)<200 %|| Cl(i,j)<0.2 || MD(i,j)>0.002 %%FA(i,j)>0.49 %MD(i,j)>0.002 % ...
                                        %( trace(i,j)<0.003 ||  ||MD(i,j)>0.002 || vol_3D(i,j)<180  )
%             QS(i,j)=0;
%          end
%     end
% end
%Use QS1 to know on which pixels Dijkstra has run 
QS1=QS;
f1=1;   % th1
f2=1;   % th2  
f3=1;   % th3
%th1 is the angle between: p(i-1),p(i),p(i+1)
th1=[];
p0=[];
tic;
imin_previous=0;
jmin_previous=0;
var1=0;
var2=0;
var3=0;
while(sum(sum(QS))>0)
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));                    
    QS(idx1(ii))=0;
    [imin,jmin]=ind2sub(size(Vx),idx1(ii));            
    %current eigenvector : 
    a_vector= [Vx(imin,jmin),Vy(imin,jmin)];
    if(sum(a_vector)>0)
        a_vector=a_vector/norm(a_vector);
    end    
    %current position index
    p1=[X(imin,jmin),Y(imin,jmin)]; 
  
    if m==1e9
        th1=[];
        var1=var1+1;
    %For the very first time that the first element doesn't have any %previous one (p(i-1))
    elseif m==0
        th1=[];   
        var2=var2+1;
    else       
        iprev=pi(imin,jmin);
        jprev=pj(imin,jmin);
        %prev position index
        %p0=[X(iprev,jprev),Y(iprev,jprev)];  
        p0=[X(iprev,jprev),Y(iprev,jprev)]-[X(imin,jmin),Y(imin,jmin)];
        var3=var3+1;
    end      
    aa=[-1,+0,+1,+1,+1,-1,+0,-1];
    bb=[+1,+1,+1,-1,+0,-1,-1,+0];
    
    for kk=1:length(aa)
        a=aa(kk);
        b=bb(kk);   
        
    %Cost calculation on the neighboring pixels 
    if (imin+a<=N) && (jmin+b<=M) && (imin+a>0) && (jmin+b>0) 
        %[X(imin+a,jmin+b),Y(imin+a,jmin+b)] is the position index of the next pixel(one out of its 8 neighbors) 
        % dp: p(i+1)-pi
        dp=( [X(imin+a,jmin+b),Y(imin+a,jmin+b)]-p1 );
        dp_nn=dp; 
        dp=dp/norm(dp);
        
        % th1
        if(~isempty(th1))
            costh1=dot( (p1-p0)/norm(p1-p0),dp/norm(dp) );
        else
            costh1=-1;
        end              
        % th2
        costh2=dot(dp,a_vector);
        % th3
        %b_vector the eigenvector on the next pixel
        b_vector=[Vx(imin+a,jmin+b),Vy(imin+a,jmin+b)];  
        if sum(b_vector)>0
            b_vector=b_vector/norm(b_vector);
        end
        costh3=dot(a_vector,b_vector);
      
        %d_1=f1*(1+costh1)+f2*(1-abs(costh2))+f3*(1-abs(costh3));
        %norm(dp_nn) euclidean distance
        %dd=(fa(imin,jmin))*d_1+(1-fa(imin,jmin))*norm(dp_nn);
        
        A=f1*(1+costh1);
        B_C=f2*(1-abs(costh2))+f3*(1-abs(costh3));
        dd= (fa(imin,jmin)*B_C)+ (1-fa(imin,jmin)*A);
        
        if d(imin+a,jmin+b)>d(imin,jmin)+dd                                         
            d(imin+a,jmin+b)=d(imin,jmin)+dd;
            pi(imin+a,jmin+b)=imin;
            pj(imin+a,jmin+b)=jmin;
        end
        imin_previous=imin;
        jmin_previous=jmin;  
    end
    end
 end
  toc  
  fprintf('\n');
end
    