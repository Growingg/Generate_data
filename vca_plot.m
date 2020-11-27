clear all;

ar = 2530;         %4991,2530
rand('state', ar);      %1996      ±¸Ñ¡ 4991  874  8013  6290  9933   7287
randn('state',234);


load 'DATA\cuprite_ref';
load 'DATA\mvcnmf_cupr_endmembers';
load 'DATA\groundTruth_Cuprite_nEnd12';
mixed = x;
E1 = M;
E1 = E1(slctBnds,:);
p = 12;
[A_est, location] = vca(mixed,'Endmembers',p,'SNR',20,'verbose','on');
li = [1,2,3,4,5,6,7,8,9,10,11,12];

figure(1)
% li = [2,5,3,4,8,6,7,1];
 for i = 1:p
    subplot(3,4,i);plot(E1(:,i),'r');hold on;plot(Aest(:,i),'g')
%     plot(A_est(:,li(i)),'b');
 end

 
figure(2)

% li = [3,8,7,2,6,11,10,1,12,5,9,4];        %4991
% li = [4,8,1,3,12,6,10,2,11,5,7,9];          %2530
% li = [1,2,3,4,5,6,7,8,9,10,11,12];          %cup2530

M1  = zeros(size(A_est));
for i =1:p
    M1(:,i) = A_est(:,li(i));
end
 for i = 1:p
    subplot(3,4,i);
    plot(E1(:,i),'r');hold on;plot(Aest(:,i),'g')
%     plot(M1(:,i),'b')
 end
 save(['.\DATA\cup\','cup_vca_endmember-2530.mat'],'M1')
