fib=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/20081006_M025Y_1Shell.src.gz.dti.mat');
trk=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/whole_brain.mat');
src=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/20081006_M025Y_1Shell.mat');
fib.dir0 = reshape(fib.dir0,[fib.dimension 3]);
fib.fa0 = reshape(fib.fa0,[fib.dimension 1]);