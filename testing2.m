%temp plotting and testing 
figure;
hold on;
[X,Y,Z]=meshgrid(1:96,1:96,1:40);
%figure;quiver3(Y,X,Z,squeeze( fib.dir0(1,:,:,:)).*fib.fa0,squeeze(fib.dir0(2,:,:,:)).*fib.fa0,squeeze(fib.dir0(1,:,:,:)).*fib.fa0);

h=plot3(Y(idx),X(idx),Z(idx),'-','Linewidth',2, 'Color',xroma);


for i=1:size(path_col,2)
    ii=round(path_lin(1,i));
    jj=round(path_col(1,i));
    kk=round(path_vol(1,i));
    
    FA=fib.fa0(ii,jj,kk);
    
    quiver3(path_lin(1,i),path_col(1,i),path_vol(1,i),...
        Vx(ii,jj,kk).*FA,...
        Vy(ii,jj,kk).*FA,...
        Vz(ii,jj,kk).*FA,4,'b');
    
end


