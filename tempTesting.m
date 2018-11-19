%plotting and tweaking
for i=1: trk.length(1)
   idxTemp(i,:)=[ trk.tracts(1,i), trk.tracts(2,i), trk.tracts(3,i) ];
   
end
figure;
plot3(idxTemp(:,1),idxTemp(:,2),idxTemp(:,3));
idxRounded=round(idxTemp);
figure;
plot3(idxRounded(:,1),idxRounded(:,2),idxRounded(:,3));title('Rounded');

figure;
plot(idxTemp(:,1),idxTemp(:,2)); title('2D');

clear i;


