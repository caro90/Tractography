%Ypologismos (gia ola ta voxels mesa sto data set)
%tou tanisti D,idiotimes,idiodianismata

%% omoia sinartisi me tn tensor_eigenVr_eigenVs_calc apla gia allou tupou dataset

%-------------------------------------------------------------------------
function [Diffusion_tensor,V,S0,Diffusion_tensor_eigenvectors,...
    Diffusion_tensor_eigenvalues,Diffusion_tensor_eigenvalues_sorted,cl,cl_normal,cp,cp_normal,cs,cs_normal,trace,FA,MD,RD]=...
    tensor_eigenVr_eigenVs_calc_version2(num_of_slices_in_diff_gradients,numOfGradients,slices_indif_gradients,gradients)

%     slices_indif_gradients=zeros(256,256,num_of_slices_in_diff_gradients,numOfGradients);
%     p=size(slices_indif_gradients,3);
%     p2=size(slices_indif_gradients,3)*size(slices_indif_gradients,4);
% 
%     temp=1;
%     for k=1:p              
%      for i=k:num_of_slices_in_diff_gradients:p2      
%             slices_indif_gradients(:,:,k,temp)=vol_3D(:,:,i);
%             temp=temp+1;
%      end
%         temp=1;
%     end

    %% xrisimopoioume ton A gia tin lisi tis eksisosis : g * D * g' = ln(S0)-ln(Sk)/b . (opou g o A) 
    A=gradients.^2; 
    A=[A,2*gradients(:,1).*gradients(:,2),2*gradients(:,1).*gradients(:,3),2*gradients(:,2).*gradients(:,3)];


    %% signal density-gia kathe slice apo to dataset me gradient=0
    S0=slices_indif_gradients(:,:,1:num_of_slices_in_diff_gradients,1);
    b=700;
    
    f1=size(slices_indif_gradients,1);
    f2=size(slices_indif_gradients,2);
    
    %% (allagi gia taxutita) o V periexei to aristero melos tis eksisosis : g * D * g' = ln (S0)-ln(Sk)/b
    V=zeros(f1,f2,num_of_slices_in_diff_gradients,numOfGradients);
        for k=1:numOfGradients
            V(:,:,:,k)=( -log( slices_indif_gradients(:,:,:,k) ) + log( S0(:,:,:) ) )/ b; 
            
        end
    
    %% o V periexei to aristero melos tis eksisosis : g * D * g' = ln (S0)-ln(Sk)/b
%     V=zeros(256,256,num_of_slices_in_diff_gradients,numOfGradients);
%     for temp=1:num_of_slices_in_diff_gradients
%         for k=1:numOfGradients
%             for i=1:256
%                for j=1:256
%            
%                 V(i,j,temp,k)=( -log( slices_indif_gradients(i,j,temp,k) ) + log( S0(i,j,temp) ) )/ b; 
%                 end
%             end
%         end
%     end

    %% Gia na mn periexei nan times to V 
            V(isnan(V))=0;
            V(isinf(V))=0;

            
    %% Ypologismos tou tanisti D

    Diffusion_tensor=zeros(f1,f2,num_of_slices_in_diff_gradients,6);

    for temp=1:num_of_slices_in_diff_gradients  % slices
        for i=1:f1
             for j=1:f2
                k=V(i,j,temp,:);
                k=reshape(k,numOfGradients,1);
                Diffusion_tensor(i,j,temp,:)=((A'*A)^-1 )*(A')* k; 
             end
        end
    end

   %% Ypologismos ton idiodianismaton kai ton idiotimon tou Diffusion tensor
    Diffusion_tensor_eigenvectors=zeros(f1,f2,num_of_slices_in_diff_gradients,9);
    Diffusion_tensor_eigenvalues=zeros(f1,f2,num_of_slices_in_diff_gradients,3);

    for y=1:num_of_slices_in_diff_gradients  %slices
      for i=1:f1
             for j=1:f2
              %temp=D_temp(i,y,:);
                temp=[Diffusion_tensor(i,j,y,1),Diffusion_tensor(i,j,y,4),Diffusion_tensor(i,j,y,5);...
                      Diffusion_tensor(i,j,y,4),Diffusion_tensor(i,j,y,2),Diffusion_tensor(i,j,y,6);...
                      Diffusion_tensor(i,j,y,5),Diffusion_tensor(i,j,y,6),Diffusion_tensor(i,j,y,3) ];
                  if ( sum(isnan(temp(:))) )>0 || (sum(isinf(temp(:)))>0 ) 
                    Diffusion_tensor_eigenvectors(i,j,y,:)=0;
            
                  else
                    [eigenvectors,a] = eig(temp);
                    a=[a(1,1),a(2,2),a(3,3)]; % o a periexei tis idiotimes     
                    Diffusion_tensor_eigenvalues(i,j,y,:)=a;   %abs(a);   ALLAGI !!!                           
                    Diffusion_tensor_eigenvectors(i,j,y,:)=eigenvectors(:);
                  end
          
             end
      end
    end

    %% Sortarima ton idiotimon

    [Diffusion_tensor_eigenvalues_sorted,index] = sort(Diffusion_tensor_eigenvalues,4,'descend');

    %% cl=linear measure , cp=planar measure , cs=spherical measure (medical image analysis 2002)
    %cl+cp+cs=1 
    % upologizo ta cp cs kai cl gia oles tis tomes 
    cl=(Diffusion_tensor_eigenvalues_sorted(:,:,:,1)-Diffusion_tensor_eigenvalues_sorted(:,:,:,2)) ./Diffusion_tensor_eigenvalues_sorted(:,:,:,1);
    cp=(Diffusion_tensor_eigenvalues_sorted(:,:,:,2)-Diffusion_tensor_eigenvalues_sorted(:,:,:,3) )./Diffusion_tensor_eigenvalues_sorted(:,:,:,1) ;
    cs=(Diffusion_tensor_eigenvalues_sorted(:,:,:,3)./Diffusion_tensor_eigenvalues_sorted(:,:,:,1));
    
    cl(isnan(cl))=0;
    cp(isnan(cp))=0;
    cs(isnan(cs))=0;
    %Allos tropos gia ton upologismo ton cp cs cl (normalized)
    
    cl_normal=(Diffusion_tensor_eigenvalues_sorted(:,:,:,1)-Diffusion_tensor_eigenvalues_sorted(:,:,:,2))./...
                sqrt(Diffusion_tensor_eigenvalues_sorted(:,:,:,1).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,2).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,3).^2);
    cp_normal=(2*(Diffusion_tensor_eigenvalues_sorted(:,:,:,2)-Diffusion_tensor_eigenvalues_sorted(:,:,:,3)))./...
                sqrt(Diffusion_tensor_eigenvalues_sorted(:,:,:,1).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,2).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,3).^2);
    cs_normal=(3*Diffusion_tensor_eigenvalues_sorted(:,:,:,3))./...
                sqrt(Diffusion_tensor_eigenvalues_sorted(:,:,:,1).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,2).^2+Diffusion_tensor_eigenvalues_sorted(:,:,:,3).^2);

    
    %% efarmogi tou index sortarismatos ton idiotimon sta idiodianismata
    for y=1:num_of_slices_in_diff_gradients %slices
        for i=1:f1
            for j=1:f2

                temp1=Diffusion_tensor_eigenvectors(i,j,y,[1: 3]);                         
                temp2=Diffusion_tensor_eigenvectors(i,j,y,[4: 6]);                         
                temp3=Diffusion_tensor_eigenvectors(i,j,y,[7: 9]);                        

                
                temp1=reshape(temp1,3,1);
                temp2=reshape(temp2,3,1);
                temp3=reshape(temp3,3,1);
                temp=[temp1,temp2,temp3];
                
                Diffusion_tensor_eigenvectors(i,j,y,1: 3 )=temp(:,index(i,j,y,1));
                Diffusion_tensor_eigenvectors(i,j,y,4: 6 )=temp(:,index(i,j,y,2));
                Diffusion_tensor_eigenvectors(i,j,y,7: 9 )=temp(:,index(i,j,y,3));
            end
        end
    end
    
    %% Trace calculation Trace(D)=?1+?2+?3
    trace=Diffusion_tensor_eigenvalues(:,:,:,1)+...
          Diffusion_tensor_eigenvalues(:,:,:,2)+...
          Diffusion_tensor_eigenvalues(:,:,:,3);     
    
   %% Fractional anisotropy 
    
    FA=sqrt(1/2)* sqrt( (Diffusion_tensor_eigenvalues(:,:,:,1)-Diffusion_tensor_eigenvalues(:,:,:,2) ).^2 + ...
                         (Diffusion_tensor_eigenvalues(:,:,:,1)-Diffusion_tensor_eigenvalues(:,:,:,3)).^2 + ...
                         (Diffusion_tensor_eigenvalues(:,:,:,2)-Diffusion_tensor_eigenvalues(:,:,:,3)).^2  ... 
                         ./(sqrt(Diffusion_tensor_eigenvalues(:,:,:,1).^2+Diffusion_tensor_eigenvalues(:,:,:,2).^2)+...
                         Diffusion_tensor_eigenvalues(:,:,:,3)).^2);
                     
  %% Mean diffusivity
    MD=( Diffusion_tensor_eigenvalues(:,:,:,1)+Diffusion_tensor_eigenvalues(:,:,:,2)+Diffusion_tensor_eigenvalues(:,:,:,3))./3;
    
  %% Radial diffusivity
    RD= (Diffusion_tensor_eigenvalues(:,:,:,2)+Diffusion_tensor_eigenvalues(:,:,:,3))./2;
    
    
end