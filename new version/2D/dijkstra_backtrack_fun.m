%%%backtrack function 

function [path_lin,path_col,B_backtrack]=dijkstra_backtrack_fun(d,pi,pj,i0,j0)
B_backtrack=zeros(size(d));
cp=1;
path_lin=[];
path_col=[];
ss=0;
while ~(pi(i0,j0)==i0 && pj(i0,j0)==j0 )
     
    i1=pi(i0,j0);
    j1=pj(i0,j0);
    i0=i1;
    j0=j1;
    ss=ss+d(i0,j0);
    path_lin(cp)=i0;
    path_col(cp)=j0;
    cp=cp+1;    
    B_backtrack(i0,j0)=1;
    if mod(cp,1000)==0
        cp
    end
end
 ss
end
