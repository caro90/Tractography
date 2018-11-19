% 3D dijkstra function version 1

function [d,p_i,pj,pz,QS]=dijkstra_function_3D(i0,j0,z0,X,Y,Z,Vx,Vy,Vz,Cl,FA,trace,MD,...
                                                                    vol_3D,SliceThickness)
global p0;
global f1;
global f2;
global f3;

%cost table:
d=ones(size(Vx))*1e9;
p_i=-ones(size(Vx));
pj=-ones(size(Vx));
pz=-ones(size(Vx));
d(i0,j0,z0)=0;
k=0;

p_i(i0,j0,z0)=i0;
pj(i0,j0,z0)=j0;
pz(i0,j0,z0)=z0;
[N,M,L]=size(p_i);
QS=ones(size(Vx));

%% trace,FA kai MD threshold gia na min ginetai dijkstra ektos kefaliou kai sto egefalonotiaio ugro
for i=1:N
    for j=1:M
        for k=1:L
            if vol_3D(i,j,k)<200  %trace(i,j,k)<0.003 || FA(i,j,k)>0.6 || MD(i,j,k)>0.002  || 
                QS(i,j,k)=0;
            end
        end
    end
end

%%
f1=1; %th1
f2=1; %th2
f3=1; %th3

%% i th1 einai i gonia metaksi ton : p(i-1),p(i),p(i+1)
tic;
th1=[];

iprev=[]; 
jprev=[];
zprev=[];
p0=[];

%%
while (sum(QS(:)>0))
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));                           
    QS(idx1(ii))=0;
    [imin,jmin,zmin]=ind2sub(size(Vx),idx1(ii));

%-----------------------------------------------------------------------------------------------  
    %current primary eigenvector : 
    a_vector= [Vx(imin,jmin,zmin),Vy(imin,jmin,zmin),Vz(imin,jmin,zmin)];
    if(norm(a_vector)>0)
        a_vector=a_vector/norm(a_vector);
    end
    %current position index
    p1=[X(imin,jmin,zmin),Y(imin,jmin,zmin),Z(imin,jmin,zmin)]; 
   
    if m==1e9      
        th1=[];     
    %gia tin proti fora pou to proto stoixeio den exei proigoumeno p(i-1)    
    elseif m==0
        th1=[];
    else
        iprev=p_i(imin,jmin,zmin);
        jprev=pj(imin,jmin,zmin);
        zprev=pz(imin,jmin,zmin);
        %prev position index:
        p0=[ X(iprev,jprev,zprev),Y(iprev,jprev,zprev),Z(iprev,jprev,zprev) ];
    end
    
    aa=[-1,+0,+1,+1,+1,-1,+0,-1,0];
    bb=[+1,+1,+1,-1,+0,-1,-1,+0,0];
    cc=[0,+1,-1];
    
    for i=1:3 % gia tin allagi ton slices
        for kk=1:length(aa)
            a=aa(kk);
            b=bb(kk);
            c=cc(i);
            
            %% katanomi ton baron(kostos) sta geitonika pixel 
            %elegxos gia na min bgoume ektos orion tis eikonas:
        if (imin+a<=N) &&  (jmin+b<=M) &&(zmin+c<=L) && (imin+a>0) &&  (jmin+b>0) && (zmin+c>0)         
            %To [X(imin+a,jmin+b),Y(imin+a,jmin+b)] einai to position index
            %tou epomenou pixel (enos apo tous 8 geitones tou current position
            %index)
        
            %I dp einai i diafora tou next me to current position index:
            %p(i+1)-pi
            dp=( [X(imin+a,jmin+b,zmin+c),Y(imin+a,jmin+b,zmin+c),Z(imin+a,jmin+b,zmin+c)]-p1 );%.*[1,1,SliceThickness];
            dp_nn=dp; 
            dp=dp/norm(dp); 
            %----------th1----------
            if(~isempty(th1))
               costh1=dot( (p1-p0)/norm(p1-p0),dp/norm(dp) ); %**** allagi gia ligoterous upologismous
                
            else
                costh1=-1;
            end
            %-----------------------
            
            %---------th2-----------
            costh2=dot(dp,a_vector);
            %-----------------------
            
            %---------th3-----------
            %to b_vector einai to idiodianisma sto epomeno pixel
            b_vector=[Vx(imin+a,jmin+b,zmin+c),Vy(imin+a,jmin+b,zmin+c),Vz(imin+a,jmin+b,zmin+c)];        
            if(norm(b_vector)>0)
                b_vector=b_vector/norm(b_vector);
            end
            costh3=dot(a_vector,b_vector);
            %------------------------
            
            d_1=f1*(1+costh1)+f2*(1-abs(costh2))+f3*(1-abs(costh3));
            %dp_nn eukleidia apostasi 
            dd=(Cl(imin,jmin,zmin))*d_1+(1-Cl(imin,jmin,zmin))*norm(dp_nn);
           
            if d(imin+a,jmin+b,zmin+c)>d(imin,jmin,zmin)+dd; 
                d(imin+a,jmin+b,zmin+c)=d(imin,jmin,zmin)+dd;
                p_i(imin+a,jmin+b,zmin+c)=imin;
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
  
  
  
    