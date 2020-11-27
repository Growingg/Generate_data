clear all
clc 
%加入程序以及数据路径
addpath(genpath('Func'));
addpath(genpath('SLIC'));
addpath(genpath('DATA'));
addpath(genpath('../results'));


%读入超像素分割的数据集
% load Indian_pines.mat
% A = indian_pines;
load jasperRidge2_R198.mat
load jasperRidge2_end4.mat
nBand = 198;
Z = reshape(Y',[nRow,nCol,nBand]);

seRadius = 1;
nIt = 10;
[row, col, D] = size(Z);

Sw = 9; Pw = round(row*col/Sw^2); Ws = 0.9;
seg = slic_HSI(Z, Pw, Ws, seRadius,nIt);
mean(seg.Cj)

nn = reshape(seg.labels,100,100);
[out] = random_color( nn,seg.labels,seg.P);
% figure(3)
% subplot_tight(1, 1, 1,[.06 .03]);
% imshow(out);

%%%%%生成丰度%%%%%
P = 4;
B = zeros(P,row*col);



Sparsity = zeros(seg.P,1);
for i=1:seg.P
    Ai = A(:,find(seg.labels==i));
    Sparsity(i) = sum(sum(Ai<0.05))/size(Ai,1)/size(Ai,2);
    Sparsity(i) = Sparsity(i)/(P/4);
%     if Sparsity(i)<0.7
%         Sparsity(i) = 0.7
%     end
    Rank = size(Ai,2);
    Bi =  Generate_abundance(P, Rank, 1, Sparsity(i));
    B(:,find(seg.labels==i))=Bi;
end

% floor(Sparsity*P)

%B为原始生成丰度, B1为丰度列之和，得到的B2为归一化之后的结果
B1 = sum(B,1);
B1 = repmat(B1,P,1);
B2 = B./B1;



% clear all
% clc 
% load Indian_pines.mat
% X = reshape(indian_pines,145*145,220);
% [mappedX, mapping] = pca(X,3);
% 
% X1 = reshape(mappedX,145,145,3);

% A =  Generate_abundance(M, P, Rank, Sparsity)