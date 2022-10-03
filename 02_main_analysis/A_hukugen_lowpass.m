%ローパスフィルタ後の主成分姿勢プロット用主成分matrix出力

%M = readmatrix('main_matrix_abs_s.csv');
M = readmatrix('output/main_matrix_lowpass.csv');
Y = readmatrix('output/main_Y.csv');
theta0 = readmatrix('output/main_theta0_abs_s.csv');

A_f = M *  Y.' + theta0;
A_fin = zeros(size(A_f,1),size(A_f,2));
for j = 1:size(A_f,2)
    A_fin(1,j) = A_f(1,j);
    for i = 2:size(A_f,1)
        A_fin(i,j) =  A_f(i,j) + A_fin(i-1,j);
    end
end
%　出力
filename = append('output/main_A_hukugen2.csv');
WRITE = A_fin.';
writematrix(WRITE,filename)