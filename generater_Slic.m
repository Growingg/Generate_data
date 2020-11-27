clear all
clc 

% ar = rand(1)*10000
ar = 6898         %3901 3178  8639   6898
rand('state', ar);      %1996      ��ѡ 4991  874  8013  6290  9933   7287
randn('state',234);

%��������Լ�����·��
addpath(genpath('Func'));
addpath(genpath('SLIC'));
addpath(genpath('DATA'));
addpath(genpath('../results'));


%���볬���طָ�����ݼ�0
% load Indian_pines.mat
% A = indian_pines;
load jasperRidge2_R198.mat
load jasperRidge2_end4.mat

% ��Ԫ����M��198*4��; ��϶�ԪY(198*10000); ��Ⱦ���A��4*10000��
nBand = 198;
Z = reshape(Y',[nRow,nCol,nBand]);

%��ԪE
P = 6;      %�ƻ����ɵļ�����Ԫ
load USGS_1995_Library.mat
wavelengths = datalib(:,1);
[dummy, indexes] = sort(wavelengths);
E1 = datalib(indexes,4:end);
names = names(4:end,:);
clear datalib;
% �����ѡP����Ԫ
indexes = randperm(size(E1,2));
indexes(12) =2 ;
E = E1(:,indexes(1:P));
[L,P] = size(E);
maxE = repmat(max(E),L,1);
minE = repmat(min(E),L,1);
E = (E-minE)./(maxE-minE);

plot(E,'DisplayName','E')

%%%���a
seRadius = 1;
nIt = 3;
[row, col, D] = size(Z);

Sw = 9; Pw = round(row*col/Sw^2); Ws = 0.9;
seg = slic_HSI(Z, Pw, Ws, seRadius,nIt);
mean(seg.Cj)

nn = reshape(seg.labels,100,100);
[out] = random_color( nn,seg.labels,seg.P);
% figure(3)
% subplot_tight(1, 1, 1,[.06 .03]);
% imshow(out);

%%%%%���ɷ��%%%%%

n = nCol*nRow;  %�ϳ����ݸ���
SNR = 30;
a = zeros(P,row*col);

%%%%%%����ÿ����Ԫ��ռ�ı���
for i=1:size(A,1)
    Prob(i) = sum(A(i,:)<0.05)/size(A,2);
end
for i = 2:size(A,1)
   Prob(i) =  Prob(i) + Prob(i-1);
end
Prob = Prob/max(Prob);

Sparsity = zeros(seg.P,1);
for i=1:seg.P
    Ai = A(:,find(seg.labels==i));
    Sparsity(i) = sum(sum(Ai~=0))/size(Ai,1)/size(Ai,2);
    if P>4
        Sparsity(i) = Sparsity(i)/(P/4);
    end
%     if Sparsity(i)<0.7
%         Sparsity(i) = 0.7
%     end
    Rank = size(Ai,2); %ÿ�����еĻ����Ԫ����
    Bi =  Generate_abundance(P, Rank, 1, Sparsity(i));
    B1(:,find(seg.labels==i))=Bi;
end



B1(B1<0.02) = 0;
% a = B1;
%BΪԭʼ���ɷ��, B1Ϊ�����֮�ͣ��õ���aΪ��һ��֮��Ľ��
B2 = sum(B1,1);
B2 = repmat(B2,P,1);
a = B1./B2;


% for i = 1:P
%     a(i,:) = smooth(a(i,:),3);
% %     A(i,:) = smooth(A(i,:));
% end


%�ϳɶ�Ԫ���Բ���
Y1 = E*a;


%�ϳɶ�Ԫ�����Բ��֣�B_aΪ˫���Ի�Ϸ�ȣ�MΪ˫���Ի�϶�Ԫ����
min_c = 0.5;
max_c = 1;
% c = min_c + (max_c - min_c)*rand(0.5*P*(P-1),n);
c = min_c;
B_a = zeros(0.5*P*(P-1),n);
step_h = 0;
for i = 1:P
    for j = i+1:P
        step_h =step_h+1;
%         c = min_c + (max_c - min_c)*rand(1,n);
        B_a(step_h,:) = a(i,:).*a(j,:);
    end
end
B = B_a.*c;

M = zeros(L,0.5*P*(P-1));
step_h = 0;
for i = 1:P
    for j = i+1:P
        step_h =step_h+1;
%         c = min_c + (max_c - min_c)*rand(1,n);
        M(:,step_h) = E(:,i).*E(:,j);
    end
end

Y2 = M*B;

%�������
y_ = (E*a).^2;
sigma = sqrt(nansum(y_(:))/n/L/10^(SNR/10));
noise = sigma*randn(L,n);
%y1Ϊ�����Ի��ģ��
y1 = Y1 + Y2;


snr =randint(L,1,[15,30]);
y = y1;
Y_a =Y1;
for i = 1:L
    y(i,:) = awgn(y1(i,:),snr(i,1));
    Y_a(i,:) = awgn(Y1(i,:),snr(i,1));
end

% y = awgn(y1,SNR);
% Y_a = awgn(Y1,SNR);
save(['.\DATA\synth\','Synthesis_data_randdB.mat'],'y')
save(['.\DATA\synth\','Synthesis_linear_data_randdB.mat'],'Y_a')
save(['.\DATA\synth\','Synthesis_data.mat'],'y1')
save(['.\DATA\synth\','Synthesis_linear_data.mat'],'Y1')
save(['.\DATA\synth\','Synthesis_linear_endmembers.mat'],'E')
save(['.\DATA\synth\','Synthesis_linear_abundance.mat'],'a')

save(['.\DATA\synth\','Synthesis_bilinear_endmembers.mat'],'M')
save(['.\DATA\synth\','Synthesis_bilinear_abundance.mat'],'B_a')
save(['.\DATA\synth\','Synthesis_bilinear_coefficient.mat'],'c')

% a1 = reshape(a,P,100,100);
% for i = 1:P
%     a11 = a1(i,:,:);
%     a11 = reshape(a11,100,100);
%     figure(i)
%     imshow(a11)
% end