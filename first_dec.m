function [k, cD, cD_L, cD_H, Sj, Sj_L, Sj_H, Sr_L, Sr_H] = first_dec(rowPpg)
%通过Tr来决定需要分解的层数，将每层细节系数分为正负，求出peak-to-sum值
%输入：原始ppg信号
%输出：论文中不难看出，留给读者自证

%------------------------------fk8小波分解----------------------------------
%小波分解
[C, L] = wavedec(rowPpg, 14, "fk8");

%第14层近似系数
cA_14 = C(1 : L(1));

%第1层到第14细节系数数组
startIndex = L(1) + 1;
for level = 2 : 1 : 15
    endPoint = startIndex + (L(level) - 1); %向量拆分的终点：起始点加数据长度
    cD{16 - level, 1} = C(startIndex : endPoint); %拆到第16 - index层时细节系数
    startIndex = endPoint + 1; %起始点迭代为终点加1处
end

%---------------------通过计算peak-to-sum ratio来确定分解层数----------------
%变量初始化
Sj = [1:14]; %每层的peak-to-sum值
Tr = 0.5; %稀疏性判定指标，可调超参数
k = 0; %分解层数

%循环遍历所有的细节系数
for level = 1 : 1 : 14
    %分解到第j层时的细节系数向量
    wj = cell2mat(cD(level)); 
    
    %分为L(负数)与H(正数)两组
    wj_L = []; %负系数数组
    wj_H = []; %正系数数组
    for index = 1 : 1 : length(wj)
        if wj(index) <= 0
            wj_L = cat(1, wj_L, wj(index));
        else
            wj_H = cat(1, wj_H, wj(index));
        end
    end

    %分别用元胞数组记录
    cD_L{level, 1} = wj_L; 
    cD_H{level, 1} = wj_H; 

    %计算peak-to-sum值，裁定第j的稀疏性
    maxwj = max(abs(wj)); %第j层最大的细节系数
    sumwj = sum(abs(wj)); %第j层细节系数绝对值之和
    Sj(level) = maxwj / sumwj;

end

%决定分解层数k
for level = 1 : 1 : 14
    if Sj(level) > Tr
        k = level;
        break
    end
end

%------------------------------计算每层Sj_L, Sj_H---------------------------


for level = 1 : 1 : 14
    maxwj_L = max(abs(cD_L{level,1})); %第j层最大的wj_L
    maxwj_H = max(abs(cD_H{level,1})); %第j层最大的wj_H
    
    sumwj_L = sum(abs(cD_L{level,1})); %第j层所有负系数绝对值相加
    sumwj_H = sum(abs(cD_H{level,1})); %第j层所有正系数绝对值相加
    
    %每层正负系数分别的peak-to-sum元胞数组
    Sj_L{level, 1} = maxwj_L / sumwj_L;
    Sj_H{level, 1} = maxwj_H / sumwj_H;

end

%------------------------------计算常值Sr_L, Sr_H---------------------------

Sr_L = (Sj_L{k, 1} + Sj_L{k+1, 1}) / 2;
Sr_H = (Sj_H{k, 1} + Sj_H{k+1, 1}) / 2;






