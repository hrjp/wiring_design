%目標姿勢計算
clear all

m = 50;         %関節数
theta_r = zeros(20,m);
filename_1 = sprintf('output/ideal_theta');

%使用する目標番号
%ideal_pat = [1 2 3 4 6 7];
%ideal_pat = [1 2 3 4 6 7 8 9 10 12 13];
%ideal_pat = [1 2 3 4 6 7 14]; 
%ideal_pat = [1 2 3 4 7 8 9 10 13];

%開きあり 0925
%ideal_pat = [1 2 3 4 7 14 8];

%1002
ideal_pat = [1 2 3 4];


pat_n = size(ideal_pat,2);      %パターン数
%theta_r 代入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1002
for k=1:1
   
    for i = 1:m
       theta_r(1,i) = 0; 
    end
    
    %wave
    %02
    %phase=[0.4,0.6,0.8,1.0];
    phase=[0.4,0.7,1.0];
    for j=2:4
        theta_r_tmp=zeros(m);
        for i=1:m
            theta_r_tmp(i)=atan2(wave_func((i)/m,phase(j-1))-wave_func((i-1)/m,phase(j-1)),1/m)*180/pi;
        end
        theta_r(j,1)=theta_r_tmp(1);
        for i=2:m
            theta_r(j,i)=theta_r_tmp(i)-theta_r_tmp(i-1);
        end
    end
    
 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%theta_r 変換
for pat = 1:pat_n
   theta_r(ideal_pat(1,pat),:) = theta_r(ideal_pat(1,pat),:) * pi/180;%[rad]
end

%出力用行列
% WRITE = zeros(m+1,pat_n);
% for pat = 1:pat_n
%     %index
%     WRITE(1,pat) = ideal_pat(1,pat);
%     for i = 1:m
%         WRITE(i+1,pat) = theta_r(ideal_pat(1,pat),i);
%     end
% end

WRITE = zeros(pat_n+1,m);
for i = 1:m
    %index
    WRITE(1,i) = i;
    for pat = 1:pat_n
        WRITE(pat+1,i) = theta_r(ideal_pat(1,pat),i);
    end
end

% WRITE = zeros(pat_n*100+1,m);
% for i = 1:m
%     %index
%     WRITE(1,i) = i;
%     for pat = 1:pat_n
%         for j = 1:100
%             WRITE((pat-1)*100+j+1,i) = theta_r(ideal_pat(1,pat),i);
%         end
%         
%     end
% end

%EXCEL出力

filename = append(filename_1,'.csv');
writematrix(WRITE,filename)
