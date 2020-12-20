%--------------------------------------------------------------------------
% This code clusters a subpicture of PaviaU datest by SPAHSIC. 
% It is for the second experiment in the paper:

% Y. Pan, Y. Jiao, T. Li and Y. Gu, "An Efficient Algorithm for 
% Hyperspectral Image Clustering", 2019 IEEE International Conference on 
% Acoustics, Speech and Signal Processing (ICASSP) 
%--------------------------------------------------------------------------

clear;close all;

%% data preprocessing
load('PaviaU.mat');
%In our experiment, we choose a subpicture of the whole PaviaU data set
HSI_3D = paviaU(151:350,101:200,:);
M = size(HSI_3D,1);
N = size(HSI_3D,2);
L = size(HSI_3D,3);
HSI = zeros(L,M*N);
% 3-D matrix to 2-D matrix
for i = 1:L
   HSI(i,:) = reshape(HSI_3D(:,:,i),1,M*N);
end
% Load the ground truth of the dataset
load('PaviaU_gt.mat');
GT = paviaU_gt;
GT = GT(151:350,101:200);
% relabel the ground truth
for i = 1:9
    if(i<3)
    [m,n] = find(GT==i);
    sta(i) = size(m,1);
    end
    if(i>3)
    m = find(GT==i);
    GT(m) = i-1;
    sta(i-1) = size(m,1);
    end
end
B_1 = GT;

%% Hyperspectral images Clustering
pre_num = 25;    % the desired number of superpixels
r = 4;           % subspace dimension
m = 0.08;        % scaling factor

[super_class, supernum] = hyperspectral_superpixels(HSI_3D,pre_num,m);        % superpixel segmentation
HSI = HSI - (mean(HSI'))'*ones(1,M*N);
U = orthvector(HSI,super_class,supernum,r);                                   % get orthonormal basis vectors of each superpixel 
aff = affinityHSI(U,supernum);                                                % get the affinity matrix
[group,C] = SpectralClustering(aff,8);                                        % spectral clustering

%% calculate performance measurements and visualization
% visualization
resultraw = zeros(1,M*N);
superraw = reshape(super_class,1,M*N);
for i = 1:supernum
    index = find(superraw == i);
    resultraw(index) = group(i);
end
result = reshape(resultraw,M,N);
result(GT==0) = 0;
figure;
imagesc(result);

% overall accuracy
GT = reshape(GT,M*N,1);
result_1=reshape(result,M*N,1);
[result_new,result_final]=bestMap(GT(GT~=0),result_1(GT~=0),result);
OA = 1-sum(result_new~=GT(GT~=0))/length(GT(GT~=0));
disp(['the overall accuracy is ', num2str(OA)])

% accuracy of each class
% acc = zeros(1,size(unique(GT(GT~=0)),1));
% for i = 1:size(unique(GT(GT~=0)),1)
%     temp = result_new(GT(GT~=0)==i);
%     num = size(find(temp==i),1);
%     acc(i) = num/size(find(GT==i),1);
% end
% acc

% kappa
for i = 1:8
    pe(i) = size(find(result_final(GT~=0)==i),1)*size(find(GT==i),1);
end
pef = sum(pe)/length(GT(GT~=0))^2;
kappa = (OA - pef)/(1-pef);
disp(['the kappa is ', num2str(kappa)])