function [new_C, L, lambda_L, lambda_H] = second_dec(k, cD, cD_L, cD_H, Sj, Sj_L, Sj_H, Sr_L, Sr_H, rowPpg1)
%求得kappa值以及阈值函数lambda并阈值去噪
%输入：第一次分解求出的各种系数数组，常量，以及ppg信号
%输出小波重组需要使用的C, L

%------------------------------fk8小波重构----------------------------------
%小波分解
[C, L] = wavedec(rowPpg1, k, "fk8");

%第k层近似系数
cA = C(1 : L(1));

%-------------------------求出各层kappa的最小值-----------------------------
kappa_L_min = [];
kappa_H_min = [];

%仅分解到第k层
for level = 1 : 1 : k
    %每一层细节系数的样本均值以及样本标准差
    miu_j = mean(cD{level, 1});
    sigma_j = std(cD{level, 1});
    
    %各层kappa的最小值
    kappa_j_L_min = (miu_j - max(cD_L{level, 1})) / sigma_j;
    kappa_j_H_min = (max(cD_H{level, 1}) - miu_j) / sigma_j;

    %填入数组
    kappa_L_min = cat(1, kappa_L_min, kappa_j_L_min);
    kappa_H_min = cat(1, kappa_H_min, kappa_j_H_min);

end


%-------------------------计算得出各层阈值-----------------------------
lambda_L = [];
lambda_H = [];

for level = 1 : 1 : k
    %每一层细节系数的样本均值以及样本标准差
    miu_j = mean(cD{level, 1});
    sigma_j = std(cD{level, 1});

    %各层kappa正负值
    kappa_j_L = ((Sr_L - Sj_H{level, 1}) / Sr_L) * kappa_L_min(level);
    kappa_j_H = ((Sr_H - Sj_H{level, 1}) / Sr_H) * kappa_H_min(level);
    
    %各层的正负阈值
    lambda_j_L = miu_j - kappa_j_L * sigma_j;
    lambda_j_H = miu_j + kappa_j_H * sigma_j;

    %填入数组
    lambda_L = cat(1, lambda_L, lambda_j_L);
    lambda_H = cat(1, lambda_H, lambda_j_H);

end

%-------------------------------阈值去噪------------------------------------
for level = 1 : 1 : k
    %对每一层的所有系数一起过阈值
    for index = 1 : 1 : length(cD{level, 1})
        if cD{level, 1}(index) >= lambda_L(level) & cD{level, 1}(index) <= lambda_H(level)
            cD{level, 1}(index) = 0;
        end
    end
end

%合成新系数向量new_C
new_C = [];
new_C = cat(1, new_C, cA);

for level = k : -1 : 1
    wj = cD{level, 1}; %分解到第j层时的细节系数向量
    new_C = cat(1, new_C, wj);
end




