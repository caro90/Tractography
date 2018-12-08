%%3D Dijkstra 
function [d,pi,pj,pz,QS]=dijkstra_function_3D_version2(i0,j0,z0,X,Y,Z,Vx,Vy,Vz,Fa)
global p0;
global f1;
global f2;
global f3;
%global f4;

%Cost table:
d=ones(size(Fa))*1e9;
pi=-ones(size(Fa));
pj=-ones(size(Fa));
pz=-ones(size(Fa));
d(i0,j0,z0)=0;
k=0;

pi(i0,j0,z0)=i0;
pj(i0,j0,z0)=j0;
pz(i0,j0,z0)=z0;
N=size(Fa,1);
M=size(Fa,2);
L=7;
QS=ones(size(Vx));

%Trace,FA kai MD thresholds to restrict the area where Dijkstra runs
% if (flag==1)    
%     for i=1:N
%         for j=1:M
%             for k=1:L
%                  if  vol_3D(i,j,k)<200 %( trace(i,j,k)<0.003 || ...
%                       %FA(i,j,k)>0.6 || MD(i,j,k)>0.002 || ...
%                       %    vol_3D(i,j,k)<200 )
%                             QS(i,j,k)=0;
%                  end
%             end
%         end
%     end
% end
f1=1;
f2=1;
f3=1;
%f4=0; 

tic;
%th1 the angle between: p(i-1),p(i),p(i+1)
th1=[];
iprev=[]; 
jprev=[];
zprev=[];
p0=[];

while (sum(QS(:)>0))
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));                           
    QS(idx1(ii))=0;
    [imin,jmin,zmin]=ind2sub(size(Vx),idx1(ii));
    
    %Current primary eigenvector
    a_vector= [Vx(imin,jmin,zmin),Vy(imin,jmin,zmin),Vz(imin,jmin,zmin)];
    
    if(norm(a_vector)>0)
        a_vector=a_vector/norm(a_vector);
    end
    %Current position index
    p1=[X(imin,jmin,zmin),Y(imin,jmin,zmin),Z(imin,jmin,zmin)]; 
   
    if m==1e9      
        th1=[];
        cos_th4_cp=[];
    %For the first time where the first element doesn't have a predecessor
    elseif m==0
        th1=[];
        cos_th4_cp=[];
    else
        iprev=pi(imin,jmin,zmin);
        jprev=pj(imin,jmin,zmin);
        zprev=pz(imin,jmin,zmin);
        %Prev position index:
        p0=[ X(iprev,jprev,zprev),Y(iprev,jprev,zprev),Z(iprev,jprev,zprev) ];
    end
    
    aa=[-1,+0,+1,+1,+1,-1,+0,-1,0];
    bb=[+1,+1,+1,-1,+0,-1,-1,+0,0];
    cc=[0,+1,-1];
    
    for i=1:3 %Slices
        for kk=1:length(aa)
            a=aa(kk);
            b=bb(kk);
            c=cc(i); 
            %Cost calculation on the neighboring voxels
        if (imin+a<=N) &&  (jmin+b<=M) &&(zmin+c<=L) && (imin+a>0) &&  (jmin+b>0) && (zmin+c>0)         
            %[X(imin+a,jmin+b),Y(imin+a,jmin+b)] is the position index
            %of the next pixel (one out of the 8 neighbors of current position index)
            %dp: p(i+1)-pi
            dp=( [X(imin+a,jmin+b,zmin+c),Y(imin+a,jmin+b,zmin+c),Z(imin+a,jmin+b,zmin+c)]-p1 );
            dp_nn=dp; 
            dp=dp/norm(dp); 
            %th1
            if(~isempty(th1))
               costh1=dot( (p1-p0)/norm(p1-p0),dp/norm(dp) );   
            else
                costh1=-1;
            end
            %th2
            costh2=dot(dp,a_vector);
            %th3
            %b_vector is the eigenvector of the next pixel/voxel
            b_vector=[Vx(imin+a,jmin+b,zmin+c),Vy(imin+a,jmin+b,zmin+c),Vz(imin+a,jmin+b,zmin+c)];        
            if(norm(b_vector)>0)
                b_vector=b_vector/norm(b_vector);
            end
            costh3=dot(a_vector,b_vector);   
            
            A=f1*(1+costh1);
            B=f2*(1-abs(costh2));
            C=f3*(1-abs(costh3));
            %dp_nn euclidean distance
            dd=Fa(imin,jmin,zmin)*(B+C)+(1-Fa(imin,jmin,zmin))*A; 
            
            if d(imin+a,jmin+b,zmin+c)>d(imin,jmin,zmin)+dd 
                d(imin+a,jmin+b,zmin+c)=d(imin,jmin,zmin)+dd;
                pi(imin+a,jmin+b,zmin+c)=imin;
                pj(imin+a,jmin+b,zmin+c)=jmin;
                pz(imin+a,jmin+b,zmin+c)=zmin;
           end
        end
      end
    end
    k=k+1;
    if mod(k,1000)==0
        fprintf('.');
    end
 end
  toc  
  fprintf('\n');
end
  
  
  
    