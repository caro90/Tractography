clear idxTemp;
% %figure;quiver3(X,Y,Z,Vx.*fib.fa0,Vy.*fib.fa0,Vz.*fib.fa0);
% figure; plot3(idxTemp(:,1),idxTemp(:,2),idxTemp(:,3)); 
% idx1=sub2ind(size(X),idxTemp(:,2),idxTemp(:,1),idxTemp(:,3));
% idx1=int32(idx1);
% hold on;% quiver3(X(idx1),Y(idx1),Z(idx1),Vx(idx1).*fib.fa0(idx1),Vy(idx1).*fib.fa0(idx1),Vz(idx1).*fib.fa0(idx1));
% quiver3(idxTemp(:,1),idxTemp(:,2),idxTemp(:,3),Vx(idx1).*fib.fa0(idx1),Vy(idx1).*fib.fa0(idx1),Vz(idx1).*fib.fa0(idx1));

figure;
fib.dir0 = reshape(fib.dir0,[3 fib.dimension]);
fib.fa0 = reshape(fib.fa0,[fib.dimension]);
Vx=squeeze(fib.dir0(1,:,:,:));
Vy=squeeze(fib.dir0(2,:,:,:));
Vz=squeeze(fib.dir0(3,:,:,:));

istart=1;
%k0 to choose the track
k0=1;
if k0>1 
    istart=sum(trk.length(1:k0-1))+1;
end
istop=istart+ trk.length(k0)-1;

%plotting from DSI datasets and tweaking
for i=istart: istop
   idxTemp(i-istart+1,:)=[ trk.tracts(1,i), trk.tracts(2,i), trk.tracts(3,i) ];
   
end

plot3(idxTemp(:,1),idxTemp(:,2),idxTemp(:,3),'Linewidth',0.5);
%idxRounded=round(idxTemp);
hold on;

for i=istart: istop
    plot3(idxTemp(i-istart+1,1),idxTemp(i-istart+1,2),idxTemp(i-istart+1,3));
    ii=round(idxTemp(i-istart+1,1));
    jj=round(idxTemp(i-istart+1,2));
    kk=round(idxTemp(i-istart+1,3));
    if kk<=0
        kk=1;
    end
    
    FA=fib.fa0(ii,jj,kk);
    
    quiver3(idxTemp(i-istart+1,1),idxTemp(i-istart+1,2),idxTemp(i-istart+1,3),...
        Vx(ii,jj,kk).*FA,...
        Vy(ii,jj,kk).*FA,...
        Vz(ii,jj,kk).*FA,4,'b');
    [FA,Vx(ii,jj,kk),Vy(ii,jj,kk),Vz(ii,jj,kk)];
end
axis equal
%plot3(idxTemp(:,2),idxTemp(:,1),idxTemp(:,3),'o-'); 