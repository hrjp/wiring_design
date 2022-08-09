%絶対座礁matrixから相対座標matrixへの変換
clear all
%主成分分析行列インポート
M = readmatrix('output/main_matrix.csv');
theta = readmatrix('output/main_theta0.csv');

M2 = zeros(size(M,1),size(M,2));
theta0 = zeros(size(theta,1),size(theta,2));

for j = 1:size(M,2)
    M2(1,j) = M(1,j);
    
    for i = 2:size(M,1)
        M2(i,j) = M(i,j)-M(i-1,j);
    end
    
end

for j = 1:size(theta,2)
    theta0(1,j) = theta(1,j);
    
    for i = 2:size(theta,1)
        theta0(i,j) = theta(i,j)-theta(i-1,j);
    end
    
end


%EXCEL出力
%data
%filename_1 = sprintf('loop_theta_r_fin_gosa_pasig_m%d_n%d_',m,n);
filename_1 = sprintf('output/main_matrix_abs_s');
filename = append(filename_1,'.csv');
writematrix(M2,filename)

filename_2 = sprintf('output/main_theta0_abs_s');
filename = append(filename_2,'.csv');
writematrix(theta0,filename)
