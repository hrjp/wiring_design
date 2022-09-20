%ローパスフィルタ後の主成分姿勢プロット
%10リンクワイヤー指モデル　準静的問題
clear all
m = 50;             %リンク数
l = zeros(1,m);     %リンク長さ [mm]
l(1,:) = 6;


l_pa = zeros(2,m);  %pa_1-p_1の距離ベクトル　l_pa(1,m)=リンク平行方向, l_pa(2,m)=リンク垂線方向 : 1からmのm個
l_pa(2,:) = 0;
l_pa0 = zeros(2,1);

l_pa(1,:) = l(1,:)/2;
l_pa0(1,1) = -1*l(1,1)/2;
l_pa0(2,1) = 0;

theta_r = zeros(1,m);
theta_h = readmatrix('output/main_A_hukugen2.csv');
theta_ideal = readmatrix('../01_ideal_pose_maker/output/ideal_theta_abs.csv');
theta_ideal(1,:) = [];
pat_n = size(theta_h,1);      %主成分要素数
color=['r','g']; %赤　理想姿勢　緑　実現姿勢
for pat = 1:pat_n
    figure
    hold on
    for state = 1:2

        file_num = sprintf('%d',pat);
        fig_file1 = append('output/hukugen_mode_',file_num,'.jpeg');
        if state==1
            theta_r = theta_ideal(pat,:);
        else
            theta_r = theta_h(pat,:);
        end

        %%
        %%plot＿理想最終状態%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        pp = zeros(2,m+2);     %リンク関節位置ベクトル 1=x, 2=y
        pr = zeros(2,m+1);     %リンク関節位置ベクトル 1=x, 2=y

        pp(1,1) = 0;
        pp(2,1) = 0;

        theta(1,1) = theta_r(1,1);

        for i = 2:m
            theta(1,i) = theta_r(1,i);
        end

            %第二リンク以降の更新
            for i = 2:m+1
                pp(1,i) = l(1,i-1) * cos(theta(1,i-1)) + pp(1,i-1);
                pp(2,i) = l(1,i-1) * sin(theta(1,i-1)) + pp(2,i-1);

                pr(1,i) = l(1,i-1) * cos(theta(1,i-1)) + pr(1,i-1);
                pr(2,i) = l(1,i-1) * sin(theta(1,i-1)) + pr(2,i-1);
            end

        x_1 = zeros(1,m+2);
        y_1 = zeros(1,m+2);

        x_1(1,1) = -l(1,1);
        y_1(1,1) = 0;

        for i=1:m+1
            x_1(1,i+1) = pp(1,i);
            y_1(1,i+1) = pp(2,i);
        end
        
        if state==2
            x_2=x_1;
            y_2=y_1;
        end
% 
% 
% 
%         tt = linspace(0,2*pi,100);
% 
%         po1 = zeros(2,m+1);
%         po2 = zeros(2,m+1);
% 
%         po1_x = zeros(1,m);
%         po1_y = zeros(1,m);
% 
%         po2_x = zeros(1,m);
%         po2_y = zeros(1,m);
% 
%         lh = 4;
% 
%         for i=1:m
% 
%             po1_x(1,i)= lh * cos(theta(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%             po1_y(1,i)= lh * sin(theta(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
% 
%             po2_x(1,i)= -1 * lh * cos(theta(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%             po2_y(1,i)= -1 * lh * sin(theta(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
% 
%             po1(1,i) = po1_x(1,i);
%             po1(2,i) = po2_x(1,i);
% 
%             po2(1,i) = po1_y(1,i);
%             po2(2,i) = po2_y(1,i);
% 
%         end
% 
%         po1(1,m+1) = -l(1,1)/2;
%         po1(2,m+1) = -l(1,1)/2;
% 
%         po2(1,m+1) = lh;
%         po2(2,m+1) = -1*lh;

        %plot(x_1,y_1,'b-o','MarkerSize',4);
        %plot(x_1,y_1,'b-','MarkerSize',1);
        plot(x_1,y_1,'g-o','LineWidth',0.5,'MarkerSize',1,'Color',color(state));
    %     for i=1:m+1
    %         hold on
    %         plot(po1(:,i),po2(:,i),'k-','LineWidth',2);
    % 
    %     end

        grid on
        %l_axis = l(1,1)*(m*2);
        l_axis = l(1,1)*(m*1.1);
        axis([-1*l_axis l_axis -1*l_axis l_axis]);
        pbaspect([1 1 1]);

        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    saveas(gcf,fig_file1)
    hold off
end
