%plotting and tweaking
for i=1: trk.length(1)
   idx(i,:)=[ trk.tracts(1,i), trk.tracts(2,i), trk.tracts(3,i) ];
   
end
figure;
plot3(idx(:,1),idx(:,2),idx(:,3));
idxNormal=round(idx);
figure;
plot3(idxNormal(:,1),idxNormal(:,2),idxNormal(:,3));

figure;
plot(idx(:,1),idx(:,2));
