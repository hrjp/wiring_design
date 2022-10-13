%相対座標の目標姿勢を絶対座標に変換
clear all

str1 = 'output/ideal_theta';
str2 = append(str1,'.csv');
A_raw = readmatrix(str2);

A_1 = A_raw(2:size(A_raw,1),:);

components = size(A_1,2);       %変数の数
number = size(A_1,1);           %データ数

A_2 = zeros(number,components);

for i = 1:number
    for j = 1:components
        
        if j ==1
             A_2(i,j) = A_1(i,j);
        else
             A_2(i,j) = A_1(i,j) + A_2(i,j-1);
        end
        
        while A_2(i,j)>pi()
            A_2(i,j) = A_2(i,j) - 2*pi();
        end
        
        while A_2(i,j)<-1*pi()
            A_2(i,j) = A_2(i,j) + 2*pi();
        end
        
    end
end

WRITE = zeros(number+1,components);
for i = 1:components
    %index
    WRITE(1,i) = i;
    for j = 1:number
        WRITE(j+1,i) = A_2(j,i);
    end
end


%EXCEL出力
filename = append(str1,'_abs.csv');
writematrix(WRITE,filename)

