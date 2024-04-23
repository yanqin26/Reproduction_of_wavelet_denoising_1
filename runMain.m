%复现论文10.1109/ACCESS.2016.2587581 作者Srivastava Madhur, Anderson C. Lindsay
%代码作者：陈羿乔 Yiqiao Chen
%日期：2024/4/23 - 2024/4/23
clc; clear; close all;
format compact; %命令行显示不换行

%==========================================================================
%-----------------------------常量设置与文件导入-----------------------------
%==========================================================================

%读取PPG文件
file1 = readmatrix("755.csv");
file1 = file1(1:length(file1),1);
rowPpg1 = cat(1, zeros(2^14 - length(file1), 1), file1); %补零延拓为2幂次的数据


%==========================================================================
%---------------------------------第一次拆分--------------------------------
%==========================================================================

%------------------------消失矩为4的多贝西平稳小波分解------------------------

%C为全尺度小波分解系数垒在一起一维系数向量
%L为记录各层系数长度的记录向量
[C, L] = wavedec(rowPpg1, 14, "db4");

%第14层近似系数
cA = C(1 : L(1));

%细节系数向量
startIndex = L(1) + 1;%向量拆分的起始点
for level = 2 : 1 : 15
    endPoint = startIndex + (L(level) - 1); %向量拆分的终点：起始点加数据长度
    cD{16 - level, 1} = C(startIndex : endPoint); %拆到第16 - index层时细节系数
    startIndex = endPoint + 1; %起始点迭代为终点加1处
end


%-------------------------通过计算peak-to-sum ratio来确定分解层数------------

%计算1到14层的细节分量稀疏性指标
old_Sj = [1 : 14];
for level = 1 : 14
    wj = cell2mat(cD(level)); %分解到第j层时的细节系数向量
    maxWj = max(abs(wj)); %最大的细节系数
    sumWj = sum(abs(wj)); %细节系数绝对值之和
    old_Sj(level) = maxWj / sumWj;
end

%弱筛选（我自己猜的啊，论文里面没有）
Tr = 0.3; %Tr: 论文中实验得出的常参是0.2，我这里目测观察决定是0.3
for level = 1 : 14
    if old_Sj(level) > Tr
        break
    end
end

k = level - 1; %k: 论文中确定分解层数的值

%绝对值最大数组
abs_max_wj_L = [];
abs_max_wj_H = [];

%计算1到k+1层的Sj_L以及Sj_H
Sj_L = [1:k+1];
Sj_H = [1:k+1];
for level = 1 : 1 : k+1
    wj = cell2mat(cD(level)); %分解到第j层时的细节系数向量
    wj_L = [];
    wj_H = [];

    %分离出正系数与负系数
    for index = 1 : 1 : length(wj)
        if wj(index) <= 0
            wj_L = cat(1, wj_L, wj(index));
        else
            wj_H = cat(1, wj_H, wj(index));
        end
    end

    %分别计算Sj
    Sj_L(level) = (max(abs(wj_L))) / sum(wj_L);
    Sj_H(level) = (max(abs(wj_H))) / sum(wj_H);

    abs_max_wj_L = cat(1, abs_max_wj_L, max(abs(Sj_L(level))));
    abs_max_wj_H = cat(1, abs_max_wj_H, max(abs(Sj_H(level))));

end

%计算常值Sr
Sr_L = (Sj_L(k) + Sj_L(k+1)) / 2;
Sr_H = (Sj_H(k) + Sj_H(k+1)) / 2;




%清除第一次拆分的变量
clear C;
clear L;
clear cD;
clear cA;

%==========================================================================
%---------------------------------第二次拆分--------------------------------
%==========================================================================

%--------------------------------拿取k层系数向量----------------------------

%C为全尺度小波分解系数垒在一起一维系数向量
%L为记录各层系数长度的记录向量
[C, L] = wavedec(rowPpg1, k, "db4");

%第k层近似系数
cA = C(1 : L(1));

%细节系数向量
startIndex = L(1) + 1;%向量拆分的起始点
for level = 2 : 1 : k+1
    endPoint = startIndex + (L(level) - 1); %向量拆分的终点：起始点加数据长度
    cD{k + 2 - level, 1} = C(startIndex : endPoint); %拆到第16 - index层时细节系数
    startIndex = endPoint + 1; %起始点迭代为终点加1处
end


%--------------------------求取各层细节系数阈值------------------------------

%计算1到k层的细节分量稀疏性指标
Sj = [1 : k];
for level = 1 : 1 : k
    wj = cell2mat(cD(level)); %分解到第j层时的细节系数向量
    maxWj = max(abs(wj)); %最大的细节系数
    sumWj = sum(abs(wj)); %细节系数绝对值之和
    Sj(level) = maxWj / sumWj;
end

%求出每层low & high 的 kappa_k_min
kappa_L_min = 1:k;
kappa_H_min = 1:k;
for level = 1 : 1 : k
    wj = cell2mat(cD(level)); %分解到第j层时的细节系数向量
    miu_j = mean(wj); %样本均值
    sigma_j = std(wj); %样本方差

    kappa_L_min(level) = (miu_j - abs_max_wj_L(level)) / sigma_j;
    kappa_H_min(level) = (abs_max_wj_H(level) - miu_j) / sigma_j;
end

%为每层进行low & high阈值去噪
for level = 1 : 1 : k
    wj = cell2mat(cD(level)); %分解到第j层时的细节系数向量
    miu_j = mean(wj); %样本均值
    sigma_j = std(wj); %样本方差
    
    %kappa值
    kappa_L = ((Sr_L - Sj_L(level)) / Sr_L) * kappa_L_min(level);
    kappa_H = ((Sr_H - Sj_H(level)) / Sr_H) * kappa_H_min(level);

    
    %阈值
    lambda_j_L = miu_j - (kappa_L * sigma_j);
    lambda_j_H = miu_j + (kappa_H * sigma_j);

    %阈值去噪
    for index = 1 : 1 : length(wj)
        if (wj >= lambda_j_L) & (wj <= lambda_j_H)
            new_cD{level, 1}(index) = 0;
        else
            new_cD{level, 1}(index) = wj(index);
        end
    end
    
end

%合成新系数向量new_C
new_C = [];
new_C = cat(1, new_C, cA);

for level = k : -1 : 1
    wj = cell2mat(new_cD(level)); %分解到第j层时的细节系数向量
    new_C = cat(1, new_C, wj');
end



%==========================================================================
%---------------------------消失矩为4的多贝西平稳小波重构--------------------
%==========================================================================

denoisiedPpg1 = waverec(new_C, L, "db4");
figure(1);
plot(rowPpg1, "b");

grid on;
hold on;

plot(denoisiedPpg1, "r");

hold off;
legend("RowPPG", "FilteredPPG")


