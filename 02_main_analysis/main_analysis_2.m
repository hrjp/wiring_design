%主成分分析プログラム
% 変数＞サンプル数に対応したver
clear all

filename_1 = sprintf('output/main_matrix');
filename_2 = sprintf('output/main_Y');
filename_3 = sprintf('output/main_avg');
filename_4 = sprintf('output/main_Y2');
filename_5 = sprintf('output/main_Y0');
filename_6 = sprintf('output/main_A_hukugen');
filename_7 = sprintf('output/main_ratios_d');
filename_8 = sprintf('output/main_theta0');

%m = 100;         %関節数
min_d = 0.7; %累積寄与率

%目標データのインポート
A_raw = readmatrix('../01_ideal_pose_maker/output/ideal_theta_abs.csv');
A = A_raw(2:size(A_raw,1),:);

components = size(A,2);       %変数の数
number = size(A,1);           %データ数

%各列要素の平均と標準偏差を求める
avg = mean(A);
s = std(A);

%読み込んだデータの各列を平均0,標準偏差1で標準化する
Normalization = zeros(number, components);      %容量の確保(計算の高速化のため)
for n = 1:number
    for comp = 1:components
        Normalization(n,comp) = ((A(n,comp) - avg(comp))/ s(comp));
    end
end

%データの相関行列を求める
R = corrcoef(A);
%相関行列の固有値の計算
[V_a,D_a] = eig(R);     % V : 固有ベクトル, D : 固有値の行列

%固有値の抽出
d_a = diag(D_a);

%並び変え dの大きい順
[d,ind] = sort(d_a,'descend');
D = D_a(ind,ind);
V = V_a(:,ind);

%寄与率
ratios_d = (d)/sum(d);

%累計寄与率
sum_ratios_d = zeros( size(d,1),1 );
sum_ratios_d(1) = ratios_d(1);
for i = 2:size(d,1)
    sum_ratios_d(i) = sum_ratios_d(i-1) + ratios_d(i);
end

%初期化
P_C_S = zeros(components,components);

% 行列の削減
j = 0;
check = 0;
for i = 1 : components    %行列Dのサイズは(componets * components)
    if check == 0       %除外条件
        j = j + 1;
        P_C_S(:,j) = V(:,i); %Principal Components Score
        if sum_ratios_d(i) > min_d  %累積寄与率がmin_dを超えるまで使う
            check = 1;
            disp(['累積寄与率 : ',num2str(sum_ratios_d(i))]);
            disp(['制御用ワイヤー本数 : ',num2str(i)]);
        end
    end
end

P_C_S = P_C_S(1:components,1:j);
Y = Normalization * P_C_S;
components2 = size(Y,2);

%確認
N2 = P_C_S *  Y.';
che_1 = Normalization - N2.';
sum_che_1 = sum(abs(che_1),2);
%

% 標準化を戻す
P_C_S_2 = zeros(components,components2);
for comp2 = 1:components2
    for comp = 1:components
        P_C_S_2(comp,comp2) = s(comp)*P_C_S(comp,comp2);
    end
end

pat_n = size(P_C_S_2,2);

%確認
A_f = P_C_S_2 *  Y.' ;
A_fin = zeros(components,number);
for i = 1:number
    A_fin(:,i) = A_f(:,i) + avg(:);
end
che_2 = A_fin - A.';
sum_che_2 = sum(abs(che_2),1);
%

%Yの負をオフセットを付ける
%成分がすべて正のY2
[Ymin,Ymin_ind] = min(Y);
Y2 = zeros(number,components2);
Ymin_2 = zeros(number,components2);
for i = 1:number
    Y2(i,:) = Y(i,:)-Ymin(1,:);
    Ymin_2(i,:) = Ymin(1,:);
end
%オフセットY0
Y0 = P_C_S_2 *  Ymin_2.';

%確認
A_f2 = P_C_S_2 *  Y2.' ;
A_fin2 = zeros(components,number);
for i = 1:number
    A_fin2(:,i) = A_f2(:,i) + Y0(:,i) + avg(:);
end
che_3 = A_fin2 - A_fin;
sum_che_3 = sum( abs(che_3) ,1);


theta0 =  zeros(components,number);

for i = 1:number
    %theta0(:,i) = Y0(:,i) + avg(:);
    theta0(:,i) = avg(:);
end



%EXCEL出力
% 変換行列 P_C_S_2
filename = append(filename_1,'.csv');
writematrix(P_C_S_2,filename)

% 主成分サンプル
filename = append(filename_2,'.csv');
WRITE = Y;
writematrix(WRITE,filename)

% A平均
filename = append(filename_3,'.csv');
WRITE = avg;
writematrix(WRITE,filename)

% Y2
filename = append(filename_4,'.csv');
WRITE = Y2;
writematrix(WRITE,filename)

% Y0
filename = append(filename_5,'.csv');
WRITE = Y0;
writematrix(WRITE,filename)

%　復元theta
filename = append(filename_6,'.csv');
WRITE = A_fin.';
writematrix(WRITE,filename)

%　寄与率　ratios_d
filename = append(filename_7,'.csv');
WRITE = ratios_d;
writematrix(WRITE,filename)

%　theta0
filename = append(filename_8,'.csv');
WRITE = theta0;
writematrix(WRITE,filename)
