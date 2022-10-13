%ワイヤリング決定

clear all

m = 30;         %関節数
l = zeros(1,m);     %リンク長さ [mm]
l(1,:) = 25;
%%n = m;             %ピン数
%pat = 5;        %目標姿勢パターン数

k = zeros(1,m);
k(1,:) = 830;
%T_max = 16.6;
T_max = 60;

%主成分分析行列インポート
M = readmatrix('../02_main_analysis/output/main_matrix_lowpass.csv');
%M = readmatrix('../02_main_analysis/output/main_matrix_abs_s.csv');
Y = readmatrix('../02_main_analysis/output/main_Y2.csv');

%目標姿勢-を排除
% Y = zeros( (size(Y_kari,1)+1)/2 ,size(Y_kari,2) );
% Y(1,:) = Y_kari(1,:);
% for i = 2:(size(Y_kari,1)+1)/2
%     Y(i,:) = Y_kari(2*(i-1),:);
% end

pat_n = size(M,2);      %主成分要素数

Y_max = max(Y,[],1);
T_abs = zeros(1,pat_n);
for i=1:pat_n
    T_abs(1,i) = T_max/Y_max(1,i);
end

%T行列の作成
T_pasig = zeros(size(Y,1),size(Y,2));
for i = 1:size(Y,1)
    for j = 1:size(Y,2)
        T_pasig(i,j) = T_abs(1,j) * Y(i,j);
    end
end


%ワイヤリング行列
A_pa = zeros(m,pat_n);
for i = 1:pat_n
    for j = 1:m
        A_pa(j,i) = M(j,i)*k(1,j)/T_abs(1,i);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% パラメータ
%トルクのつり合い式
%m = 100;             %リンク数
% k = zeros(1,m);
% k(1,:) = 100;

l_pa = zeros(2,m);  %pa_1-p_1の距離ベクトル　l_pa(1,m)=リンク平行方向, l_pa(2,m)=リンク垂線方向 : 1からmのm個
l_pa(2,:) = 0;
l_pa0 = zeros(2,pat_n);

l_pa(1,:) = l(1,:)/2;
l_pa0(1,:) = -1*l(1,1)/2;
l_pa0(2,:) = 0;
% l_pa0(2,1) = -0.5;
% l_pa0(2,2) = 1;

% pasig計算変数
p = zeros(2,m+1);     %リンク関節位置ベクトル 1=x, 2=y
theta = zeros(1,m); %リンク関節角度[rad]
thetasig = zeros(1,m); %リンク関節角度[rad]
pa = zeros(2,m+1);    %リンクワイヤー通し中心ベクトル 1=x, 2=y : 0からmのm+１個
eb = zeros(2,m+1);     %張力tb単位ベクトル : 0からmのm+１個
K = zeros(m,m);     %ばねトルク係数項
pasig = zeros(2,m);

% theta_r = zeros(2,m);
% for 
%     theta_r()
% end


% エラー記録
ernum = 1;
error_pin = zeros(ernum,m);
error_num =0;
pin_error = 0;      %pinの解の有無


% 記録用
pasig_loop = zeros(2,m,pat_n);
WRITE = zeros(3*pat_n,m+1);
plot_pasig_1 = zeros(2,m+1);
plot_pasig_2 = zeros(2,m+1);

%%
% pasig 計算
for loop = 1:pat_n
    % n_lp番ワイヤーについて
    pa0 = zeros(2,1);
    pa0(1,1) = l_pa0(1,loop);
    
    l_pa0(2,loop) = A_pa(1,loop);
    pa0(2,1) = l_pa0(2,loop);
    
    pasig(1,1) = l_pa(1,1);

    % theta = 0 とみなしたときの近似
    pa0sig(1,1) = pa0(1,1);
    pa0sig(2,1) = pa0(2,1);

    A_1 = ( A_pa(1,loop) )^2 - pa0sig(1,1)^2;
    A_2 = -2*( (A_pa(1,loop))^2 )*pa0sig(2,1) + 2*pa0sig(1,1)*pa0sig(2,1)*pasig(1,1); 
    A_3 = ( (A_pa(1,loop))^2 ) * ( (pa0sig(1,1)-pasig(1,1))^2 + pa0sig(2,1)^2 ) - (pasig(1,1)^2) * (pa0sig(2,1)^2);

    pasig_check = A_2^2 - 4*A_1*A_3;

    %確認用　解の存在条件
    Z = A_pa(1,loop);
    X1 = pasig(1,1);
    X0 = pa0sig(1,1);
    Y0 = pa0sig(2,1);

    D = 4*Z^2 * (X0-X1)^2 * ( -1*Z^2 + ( X0^2 + Y0^2 ) );

    Z_hani = ( X0^2 +Y0^2 )^(1/2);
    
    pasig_k1 = (-1 * A_2 + sqrt( A_2^2 - 4*A_1*A_3 ))/(2*A_1);
    pasig_k2 = (-1 * A_2 - sqrt( A_2^2 - 4*A_1*A_3 ))/(2*A_1);

%     if abs(pasig_k1) < 10^(-5) 
%         pasig_k1 = 0;
%     end 
% 
%     if abs(pasig_k2) < 10^(-5) 
%         pasig_k2 = 0;
%     end 
    %pasig　正負チェック
    pasig_e_1 = pasig(1,1)*(pa0sig(2,1)-pasig_k1) - pasig_k1*(pa0sig(1,1)-pasig(1,1));
    pasig_e_2 = pasig(1,1)*(pa0sig(2,1)-pasig_k2) - pasig_k2*(pa0sig(1,1)-pasig(1,1));
    
    if isreal(pasig_k1) == 0
        %pasig_k1が虚数
        if isreal(pasig_k2) == 0
            %pasig_k1，pasig_k2は虚数
            fprintf('解なし 1_1\n');

            pin_error = 1;
            error_num = error_num + 1;
            for er = 1:m
                error_pin(error_num,er) = er;
            end
            %return
        else
             %pasig_k1は虚数，pasig_k2が実数
            if pasig_e_2 * A_pa(1,loop) <0
                %pasig_k2とthetaが同方向ではない
                fprintf('解なし 1_2\n');
                %pn
                pin_error = 1;
                error_num = error_num + 1;
                for er = 1:m
                    error_pin(error_num,er) = er;
                end
                %return
            else
                %解あり　pasig_k2
                pasig(2,1) = pasig_k2;
            end
        end
    else
        %pasig_k1は実数
        if isreal(pasig_k2) == 0
            %pasig_k1は実数，pasig_k2は虚数
            if pasig_e_1 * A_pa(1,loop) < 0
                %pasig_k1とthetaが同方向でない
                fprintf('解なし 1_3\n');
                %pn
                pin_error = 1;
                error_num = error_num + 1;
                for er = 1:m
                    error_pin(error_num,er) = er;
                end
                %return
            else
                pasig(2,1) = pasig_k1;
            end
        else
            %pasig_k1，pasig_k2が実数
            if pasig_e_1 * A_pa(1,loop) < 0 && pasig_e_2 * A_pa(1,loop) < 0
                %pasig_k1，pasig_k2とthetaが同方向でない
                fprintf('解なし 1_4\n');
                %pn
                pin_error = 1;
                error_num = error_num + 1;
                for er = 1:m
                    error_pin(error_num,er) = er;
                end
                %return
            else
                %pasig_k1もしくはpasig_k2がthetaと同方向
%                 if pasig_e_1 * A_pa(1,loop) < 0 && pasig_e_2 * A_pa(1,loop) < 0
%                     %pasig_k1とpasig_k2がthetaと同方向でない
%                     fprintf('解なし 1_5\n');
%                     %pn
%                     pin_error = 1;
%                     error_num = error_num + 1;
%                     for er = 1:m
%                         error_pin(error_num,er) = er;
%                     end
%                     %return
%                 else
                    if pasig_e_1 * A_pa(1,loop) >= 0 && pasig_e_2 * A_pa(1,loop) >= 0
                        %pasig_k1とpasig_k2がthetaと同方向
                        if abs(pasig_k1 - pa0sig(2,1)) > abs(pasig_k2 - pa0sig(2,1))
                            pasig(2,1) = pasig_k2;
                        else
                            pasig(2,1) = pasig_k1;
                        end
                    else
                        %pasig_k1もしくはpasig_k2がthetaと同方向
                        if pasig_e_1 * A_pa(1,loop) >= 0
                            pasig(2,1) = pasig_k1;
                        else
                            pasig(2,1) = pasig_k2;
                        end
                    end
                %end
            end
        end
    end


    for ip = 2:m
        pasig(1,ip) = l_pa(1,ip);

        % theta = 0 とみなしたときの近似
        pa0sig(1,1) = pasig(1,ip-1) - l(1,ip-1);
        pa0sig(2,1) = pasig(2,ip-1);
        
        A_1 = ( A_pa(ip,loop) )^2 - pa0sig(1,1)^2;
        A_2 = -1*( A_pa(ip,loop) )^2 * 2 * pa0sig(2,1) + 2 * pa0sig(1,1) * pa0sig(2,1) * pasig(1,ip); 
        A_3 = ( A_pa(ip,loop) )^2 * ( (pa0sig(1,1)-pasig(1,ip))^2 + pa0sig(2,1)^2 ) - pasig(1,ip)^2 * pa0sig(2,1)^2;

        pasig_k1 = (-1 * A_2 + sqrt( A_2^2 - 4*A_1*A_3 ))/(2*A_1);
        pasig_k2 = (-1 * A_2 - sqrt( A_2^2 - 4*A_1*A_3 ))/(2*A_1);

%         if abs(pasig_k1) < 10^(-5) 
%             pasig_k1 = 0;
%         end 
% 
%         if abs(pasig_k2) < 10^(-5) 
%             pasig_k2 = 0;
%         end 
        %pasig　正負チェック
        pasig_e_1 = pasig(1,ip)*(pa0sig(2,1)-pasig_k1) - pasig_k1*(pa0sig(1,1)-pasig(1,ip));
        pasig_e_2 = pasig(1,ip)*(pa0sig(2,1)-pasig_k2) - pasig_k2*(pa0sig(1,1)-pasig(1,ip));
    
        
        %pasig(2,i) = pasig_k2; 
        if isreal(pasig_k1) == 0
            %pasig_k1が虚数
            if isreal(pasig_k2) == 0
                %pasig_k1，pasig_k2は虚数
                fprintf('解なし %d_1\n',ip);
                %pn
                pin_error = 1;
                error_num = error_num + 1;
                for er = 1:m
                    error_pin(error_num,er) = er;
                end
                %return
            else
                 %pasig_k1は虚数，pasig_k2が実数
                if pasig_e_2 * A_pa(ip,loop) <0
                    %pasig_k2とthetaが同方向ではない
                    fprintf('解なし %d_2\n',ip);
                    %pn
                    pin_error = 1;
                    error_num = error_num + 1;
                    for er = 1:m
                        error_pin(error_num,er) = er;
                    end
                    %return
                else
                    %解あり　pasig_k2
                    pasig(2,ip) = pasig_k2;
                end
            end
        else
            %pasig_k1は実数
            if isreal(pasig_k2) == 0
                %pasig_k1は実数，pasig_k2は虚数
                if pasig_e_1 * A_pa(ip,loop) < 0
                    %pasig_k1とthetaが同方向でない
                    fprintf('解なし %d_3\n',ip);
                    %pn
                    pin_error = 1;
                    error_num = error_num + 1;
                    for er = 1:m
                        error_pin(error_num,er) = er;
                    end
                    %return
                else
                    pasig(2,ip) = pasig_k1;
                end
            else
                %pasig_k1，pasig_k2が実数
                if pasig_e_1 * A_pa(ip,loop) < 0 && pasig_e_2 * A_pa(ip,loop) < 0
                    %pasig_k1，pasig_k2とthetaが同方向でない
                    fprintf('解なし %d_4\n',ip);
                    %pn
                    pin_error = 1;
                    error_num = error_num + 1;
                    for er = 1:m
                        error_pin(error_num,er) = er;
                    end
                    %return
                else
                    %pasig_k1もしくはpasig_k2がthetaと同方向
                    if pasig_e_1 * A_pa(ip,loop) >= 0 && pasig_e_2 * A_pa(ip,loop) >= 0
                        %pasig_k1とpasig_k2がthetaと同方向
                        if abs(pasig_k1 - pa0sig(2,1)) > abs(pasig_k2 - pa0sig(2,1))
                            pasig(2,ip) = pasig_k2;
                        else
                            pasig(2,ip) = pasig_k1;
                        end
                    else
                        %pasig_k1もしくはpasig_k2がthetaと同方向
                        if pasig_e_1 * A_pa(ip,loop) >= 0
                            pasig(2,ip) = pasig_k1;
                        else
                            pasig(2,ip) = pasig_k2;
                        end
                    end
                    
                end
            end
        end
    end
    %pasig 記録
    pasig_loop(:,:,loop) = pasig(:,:);
    
    %書き込み
    WRITE( (3*loop-2) ,1) = loop;
    WRITE( (3*loop-1) ,1) = l_pa0(1,loop);
    WRITE( (3*loop) ,1) = l_pa0(2,loop);
    WRITE( (3*loop-1):(3*loop) ,2:m+1) = pasig_loop(:,:,loop);
end

plot_pasig_1(1,1) = 0;
plot_pasig_1(2,1) = l_pa0(2,1);
plot_pasig_2(1,1) = 0;
plot_pasig_2(2,1) = l_pa0(2,2);
    
for i = 2:m+1
    plot_pasig_1(1,i) = i-1;
    plot_pasig_1(2,i) = pasig_loop(2,i-1,1);
    plot_pasig_2(1,i) = i-1;
    plot_pasig_2(2,i) = pasig_loop(2,i-1,2);
end


figure
plot(plot_pasig_1(1,:),plot_pasig_1(2,:),'r');
hold on
plot(plot_pasig_2(1,:),plot_pasig_2(2,:),'b');

%axis([ 0 50 -4 7]);

%EXCEL出力
%data
%filename_1 = sprintf('loop_theta_r_fin_gosa_pasig_m%d_n%d_',m,n);
filename_1 = sprintf('output/pasig_loop2');
filename = append(filename_1,'.csv');
writematrix(WRITE,filename)

filename_2 = sprintf('output/pasig_T');
filename = append(filename_2,'.csv');
writematrix(T_pasig,filename)
