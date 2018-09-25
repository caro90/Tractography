%%%2D dijkstra function

function [d,pi,pj,QS]=dijkstra_testing_with_VectorsFields_djk_fun(i0,j0,X,Y,Vx,Vy,Cl)
global p0;
global f1;
global f2;
global f3;
global f4;

%[N,M]=size(Vx);

%cost table:
d=ones(size(Vx))*1e9;
pi=-ones(size(Vx));
pj=-ones(size(Vx));
d(i0,j0)=0;
k=0;

pi(i0,j0)=i0;
pj(i0,j0)=j0;

QS=ones(size(Vx));
%m0=max(max(d));
f1=0.1;
f2=0.1;
f3=1;
f4=1;
tic;
while(sum(sum(QS))>0)
    
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));                    
    QS(idx1(ii))=0;
    [imin,jmin]=ind2sub(size(Vx),idx1(ii));
%-----------------------------------------------------------------------------------------------  
    %f=atan2(Y(imin,jmin),X(imin,jmin));             

    %current eigenvector : 
    a_vector= [Vx(imin,jmin),Vy(imin,jmin)];
    if(sum(a_vector)>0)
        a_vector=a_vector/norm(a_vector);
    end
    %current position index
    p1=[X(imin,jmin),Y(imin,jmin)]; 
   
    %i th1 einai i gonia metaksi ton : p(i-1),p(i),p(i+1)
    th1=[];
    iprev=[];
    jprev=[];
    p0=[];
    if m==1e9  %???pote einai to m=1e9?axreiasto if
        th1=[];
       
    %gia tin proti fora pou to proto stoixeio den exei proigoumeno p(i-1)    
    elseif m==0
        th1=[];
       
    else
        iprev=pi(imin,jmin);
        jprev=pj(imin,jmin);
        %prev position index :
        p0=[X(iprev,jprev),Y(iprev,jprev)];        
       
    end
    
    aa=[-1,+0,+1,+1,+1,-1,+0,-1];
    bb=[+1,+1,+1,-1,+0,-1,-1,+0];
    
    for kk=1:length(aa)
        [d,pi,pj]=my_fun(imin,jmin,p1,th1,X,Y,Vx,Vy,d,pi,pj,aa(kk),bb(kk),a_vector,Cl);
    end
  
    k=k+1;
    if mod(k,1000)==0
        fprintf('.');
    end
end
  toc  
  fprintf('\n');
end

%% katanomi ton baron(kostos) sta geitonika pixel 
function [d,pi,pj]=my_fun(imin,jmin,p1,th1,X,Y,Vx,Vy,d,pi,pj,a,b,a_vector,Cl)
global p0;
global f1;
global f2;
global f3;
global f4;

[N,M]=size(pi);
    %elegxos gia na min bgoume ektos orion tis eikonas:
    if (imin+a<=N) &&  (jmin+b<=M) && (imin+a>0) &&  (jmin+b>0) 
        %To [X(imin+a,jmin+b),Y(imin+a,jmin+b)] einai to position index
        %tou epomenou pixel (enos apo tous 8 geitones tou current position
        %index)
        
        %I dp einai i diafora tou next me to current position index:
        %p(i+1)-pi
        dp=( [X(imin+a,jmin+b),Y(imin+a,jmin+b)]-p1 );
        dp_nn=dp; 
        dp=dp/norm(dp);  
        if(~isempty(th1))
            costh1=dot( (p1-p0)/norm(p1-p0),dp/norm(dp) );
        else
            costh1=-1;
        end                
        costh2=dot(dp,a_vector);
        %to b_vector einai to idiodianisma sto epomeno pixel
        b_vector=[Vx(imin+a,jmin+b),Vy(imin+a,jmin+b)];  
        if sum(b_vector)>0
            b_vector=b_vector/norm(b_vector);
        end
        costh3=dot(a_vector,b_vector);
                              
        dd=f1*(1+costh1)+f2*(1-abs(costh2))+f3*(1-abs(costh3));
        
        %dp_nn eukleidia apostasi 
        
        dd=Cl(imin,jmin)*dd+(1-Cl(imin,jmin))*norm(dp_nn)*f4;
        
        if d(imin+a,jmin+b)>d(imin,jmin)+dd; 
            d(imin+a,jmin+b)=d(imin,jmin)+dd;
            pi(imin+a,jmin+b)=imin;
            pj(imin+a,jmin+b)=jmin;
        end
    end
end
    