%%%2D dijkstra function

%version 11/2/17

function [d,p_i,pj,QS1]=dijkstra_function(i0,j0,X,Y,Vx,Vy,Cl,trace,FA,MD,vol_3D,cp,flag)
global p0;
global f1;
global f2;
global f3;
global f4;

%% cost table:
d=ones(size(Vx))*1e9;
d(i0,j0)=0;

%% Stous pinakes p_i kai pj (previous) apothikeuoume kathe fora ton proigoumeno pixel komvo (position index) 
p_i=-ones(size(Vx));
pj=-ones(size(Vx));
%k=0;
%Thetoume tis sintetagmenes tou arxikou simeiou stous pinakes p_i kai pj
%gia na to xrisimopoihsoume gia seed sto backtrack :
p_i(i0,j0)=i0;
pj(i0,j0)=j0;

%% Sto pinaka Qs thetoume 0 sto pixel-komvo pou epeksergastikame 
QS=ones(size(Vx));
[N,M]=size(p_i);
%% trace,FA kai MD threshold gia na min ginetai dijkstra ektos kefaliou kai sto egefalonotiaio ugro

if flag==1
    for i=1:N
        for j=1:M
            if vol_3D(i,j)<200 || trace(i,j)>0.008
                %|| Cl(i,j)<0.2 || MD(i,j)>0.002 ...
                %  FA(i,j)>0.49 || MD(i,j)>0.002 ...
                % trace(i,j)<0.003 || MD(i,j)>0.002 ...
                %|| vol_3D(i,j)<180  
                QS(i,j)=0;
             end
        end
    end
end
QS1=QS;
%%
f1=1;   % th1
f2=1;   % th2  
f3=1;   % th3
f4=5;

%% i th1 einai i gonia metaksi ton : p(i-1),p(i),p(i+1)
th1=[];
p0=[];

tic;
imin_previous=0;
jmin_previous=0;
while(sum(sum(QS))>0)
    
    idx1=find(QS==1);
    [m,ii]=min(d(idx1));                    
    QS(idx1(ii))=0;
    [imin,jmin]=ind2sub(size(Vx),idx1(ii));
           
    %current eigenvector : 
    a_vector= [Vx(imin,jmin),Vy(imin,jmin)];
    %Gia na apotrepetai i diairesi me to 0
    if(norm(a_vector)>0)
        a_vector=a_vector/norm(a_vector);
    end
        
    %current position index
    p1=[X(imin,jmin),Y(imin,jmin)]; 
   
    if m==1e9
        th1=[];
    %gia tin proti fora pou to proto 
    %stoixeio den exei proigoumeno p(i-1)
    elseif m==0
        th1=[];   
    else       
        th1=0;
        iprev=p_i(imin,jmin);
        jprev=pj(imin,jmin);
        %prev position index : [X(iprev,jprev),Y(iprev,jprev)]
        %p0 i diafora tou previous me to current position index p(i-1) - pi :
        p0= [X(iprev,jprev),Y(iprev,jprev)]-[X(imin,jmin),Y(imin,jmin)] ;
        p0=p0/norm(p0);
    end
      
    aa=[-1,+0,+1,+1,+1,-1,+0,-1];
    bb=[+1,+1,+1,-1,+0,-1,-1,+0];
    
    t=zeros(1,9);
    for kk=1:length(aa)
        a=aa(kk);
        b=bb(kk);   
        
    %% katanomi ton baron(kostos) sta geitonika pixel 
    %elegxos gia na min bgoume ektos orion tis eikonas:
    if (imin+a<=N) &&  (jmin+b<=M) && (imin+a>0) &&  (jmin+b>0) 
        
        %%
        %To [X(imin+a,jmin+b),Y(imin+a,jmin+b)] einai to position index
        %tou epomenou pixel (enos apo tous 8 geitones tou current position
        %index)
        
        % I dp einai i diafora tou next 
        % me to current position index:
        % p(i+1)-pi
        dp=( [X(imin+a,jmin+b),Y(imin+a,jmin+b)]-p1 );
        dp_nn=dp; 
        dp=dp/norm(dp);
        %% ------th1-------
        if(~isempty(th1))
            costh1=dot( p0,dp)/norm(p0)*norm(dp);
            
        else
            costh1=-1;
        end  
        %----------------------
        
        %% --------th2----------
        costh2=dot(dp,a_vector);
        %----------------------
        
        %% -----------th3--------
        %to b_vector einai to idiodianisma sto epomeno pixel
        b_vector=[Vx(imin+a,jmin+b),Vy(imin+a,jmin+b)];  
        
        %Gia na apotrepetai i diairesi me to 0
        if norm(b_vector)>0 
            b_vector=b_vector/norm(b_vector);
        end
        costh3=dot(a_vector,b_vector);
        %-------------------------
        
        d_1=f1*(1+( costh1 ))+f2*(1-abs(costh2))+f3*(1-abs(costh3));
       
        %norm(dp_nn) eukleidia apostasi 
        %dd=(Cl(imin,jmin))*d_1 + f4*(1-Cl(imin,jmin)).^2*norm(dp_nn);

        dd=(1-(1-Cl(imin,jmin)).^3)*d_1 + f4*(1-Cl(imin,jmin)).^3*norm(dp_nn);
        
        %% *** 
        %dd=(Cl(imin,jmin))*d_1  +norm(dp_nn)* Cp(imin,jmin);
%        if(L1(imin,jmin)==0)
%             dd=1000;
%        else
%             dd=dd*L1max/(L1(imin,jmin));
%        end

      
        if d(imin+a,jmin+b)>d(imin,jmin)+dd 
            d(imin+a,jmin+b)=d(imin,jmin)+dd;
            p_i(imin+a,jmin+b)=imin;
            pj(imin+a,jmin+b)=jmin;
        end
        imin_previous=imin;
        jmin_previous=jmin;
        
    end
    end

 end
  toc   
end
    