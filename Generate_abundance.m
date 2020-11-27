function A = Generate_abundance(M, P, Rank, Sparsity)
if Rank < 0 || Rank > P || Sparsity < 0 || Sparsity >= 1
    error('The Rank or Sparsity is wrong');
end
% Rank多少种混合模式相同，Sparsity每次有（1-Sparsity）*10个端元混合  P为生成混合像元的个数 M为由M个端元混合
Sparsity = 1-Sparsity;
A = rand(M,P);
Res = mod(P, Rank);
Num = floor(P/Rank);
Sn = round(Sparsity*M);
if Sn==M-1
    Sn=M-2;
end
if Res ~= 0
    for i=1:Rank-1
        I = randperm(M,Sn);
        tt = rand(M,1).*rand(M,1);
        zz = zeros(M,Num);
        for j = 1:M
            zz(j,:) = normrnd(tt(j),0.03,1,Num);
%             zz(j,:) =tt(j)+ 0.05*rand(1,Num);
        end
        A(:,1+(i-1)*Num:i*Num) = zz;
%         A(:,1+(i-1)*Num:i*Num) = rand(M,Num).*rand(M,Num);
        A(I,1+(i-1)*Num:i*Num) = 0;
    end
    I = randperm(M,Sn);
    tt = rand(M,1).*rand(M,1);
    zz = zeros(M,Num + Res);
    for j = 1:M
        zz(j,:) = normrnd(tt(j),0.03,1,Num + Res);
%         zz(j,:) =tt(j)+ 0.05*rand(1,Num + Res);
    end
     A(:,1+(Rank-1)*Num:end) = zz;
%     A(:,1+(Rank-1)*Num:end) = rand(M,Num + Res).*rand(M,Num + Res);
    A(I,1+(Rank-1)*Num:end) = 0;
else
    for i=1:Rank
        I = randperm(M,Sn);
        tt = mean(rand(M,1).*rand(M,1),2);
        zz = zeros(M,Num);
        for j = 1:M
            zz(j,:) = normrnd(tt(j),0.03,1,Num);
%             zz(j,:) =tt(j)+ 0.05*rand(1,Num);
        end
        A(:,1+(i-1)*Num:i*Num) = zz;
%         A(:,1+(i-1)*Num:i*Num) = rand(M,Num).*rand(M,Num);
        A(I,1+(i-1)*Num:i*Num) = 0;
    end
     
end
A = abs(A);
A1 = sum(A,1);
A2 = repmat(A1,M,1);
A = A./A2;

