%目標姿勢計算
clear all

m = 50;         %関節数
theta_r = zeros(20,m);
filename_1 = sprintf('ideal_theta_1026_1');

%使用する目標番号
%ideal_pat = [1 2 3 4 6 7];
%ideal_pat = [1 2 3 4 6 7 8 9 10 12 13];
%ideal_pat = [1 2 3 4 6 7 14]; 
%ideal_pat = [1 2 3 4 7 8 9 10 13];

%開きあり 0925
%ideal_pat = [1 2 3 4 7 14 8];

%1002
ideal_pat = [1 2 3 4 5];


pat_n = size(ideal_pat,2);      %パターン数
%theta_r 代入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1002
for k=1:1
   
    for i = 1:m
       theta_r(1,i) = 0; 
    end

    %02
    i2 = 60/100*m;
    for i = 1:i2
       theta_r(2,i) = 1.2/i2*i*100/m; 
    end
    for i = i2+1:m
       theta_r(2,i) = 1.6/(m-i2)*i*100/m; 
    end
    
    %03 
    i2 = 20/100*m;
    for i = 1:i2
       theta_r(3,i) = -40/(i2); 
    end
    for i = i2+1:m
       theta_r(3,i) = 160/(m-i2); 
    end    
    
    %04
    i1 = round(m/6);
    i2 = round(i1*2);
    i3 = round(i1*3.5);
    for i = 1:i2
       theta_r(4,i) = 40/i2; 
    end
    for i = i2+1:i3
       theta_r(4,i) = -100/(i3-i2);
    end
    for i = i3+1:m
       theta_r(4,i) = 135/(m-i3);
    end
        
    %05 開き
    for i = 1:m
       theta_r(5,i) = -2/m*i*100/m;  %01
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
