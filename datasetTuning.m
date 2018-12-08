fib=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/20081006_M025Y_1Shell.src.gz.dti.mat');
trk=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/whole_brain.mat');
src=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/20081006_M025Y_1Shell.mat');
fib.dir0 = reshape(fib.dir0,[3 fib.dimension]);
fib.fa0 = reshape(fib.fa0,[fib.dimension]);
%Old Dataset Larissa's Hospital
oldData=load('C:\Users\Chris\Documents\MATLAB\Tractography\data/DTI_dataset_egefalografias_Nosokomeio_Larissas.mat');