%%Backtrack 3D

function [path_lin,path_col,path_vol,B_backtrack]=dijkstra_backtrack_fun_3D(d,pi,pj,pz,i0,j0,z0)
    B_backtrack=zeros(size(d));
    cp=1;
    path_lin=[];
    path_col=[];
    path_vol=[];
    ss=0;
    path_lin(cp)=i0;
    path_col(cp)=j0;
    path_vol(cp)=z0;
    cp=cp+1; 
    while ~(pi(i0,j0,z0)==i0 && pj(i0,j0,z0)==j0 && pz(i0,j0,z0)==z0)
        i1=pi(i0,j0,z0);
        j1=pj(i0,j0,z0);
        z1=pz(i0,j0,z0);
        i0=i1;
        j0=j1;
        z0=z1;
        ss=ss+d(i0,j0,z0);
        path_lin(cp)=i0;
        path_col(cp)=j0;
        path_vol(cp)=z0;
        cp=cp+1;    
        B_backtrack(i0,j0,z0)=255;
        if mod(cp,1000)==0
            cp
        end
    end
end