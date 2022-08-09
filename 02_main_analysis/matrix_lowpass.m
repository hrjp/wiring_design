%ローパスフィルタ

M = readmatrix('output/main_matrix_abs_s.csv');
M2 = zeros(size(M,1),size(M,2));
x = zeros(size(M,1),size(M,2));

for i = 1:size(M,2)
    x(:,1) = M(:,i);

    %lowpassフィルタ
%     ws = 0.25;
%     M2(:,i) = lowpass(x(:,1),ws);
    
    %移動平均
%     windowSize = 5; 
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     M2(:,i) = filter(b,a,x(:,1));
    
    %バタワース
    n=2;
    Wn=0.2;
    x2 = zeros(size(M,1),1);
    for j= 1:size(M,1)
        x2(j,1) = x(j,1)-x(1,1);
    end
    [b,a] = butter(n,Wn);
    x3 = filter(b,a,x2(:,1));
    x4 = flip(x3);
    x5 = zeros(size(M,1),1);
    for j= 1:size(M,1)
        x5(j,1) = x4(j,1)-x4(1,1);
    end
    x6 = filter(b,a,x5);
    for j= 1:size(M,1)
        x6(j,1) = x6(j,1) + x4(1,1);
    end
    
    M2(:,i) = flip(x6);
    
    for j= 1:size(M,1)
        M2(j,i) = M2(j,i) + x(1,1);
    end
    
    
    figure
    plot(M(:,i),'g')
    hold on 
    plot(M2(:,i),'r')
    hold off
end



filename = sprintf('output/main_matrix_lowpass.csv');
WRITE = M2;
writematrix(WRITE,filename)