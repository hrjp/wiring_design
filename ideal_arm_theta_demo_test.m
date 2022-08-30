clear all
% 
% num = 4;
% 
% %m = 100;
% m = 50;             %link size
% l = zeros(1,m);     
% %l(1,:) = 2;
% l(1,:) = 6;
% lh = 10;     %link width
% %k = zeros(1,m);
% %k(1,:) = 100;
% l_pa = zeros(2,m);
% l_pa(2,:) = 0;
% l_pa0 = zeros(2,1);
% 
% l_pa(1,:) = l(1,:)/2;
% l_pa0(1,1) = -1*l(1,1)/2;
% l_pa0(2,1) = 0;
% 
% %str1 ="kaiseki_m50_1deg_miss_0831";
% %fig_file2 = append('ideal_mode_0921_07.jpeg');
% fig_file2 = append('output/ideal_',string(num),'.jpeg');
% theta_r = zeros(20,m);
% 
% 
%     %01
%     for i = 1:m
%        theta_r(1,i) = 0; 
%     end
%     
%    %wave
%     %02
%     %phase=[0.4,0.6,0.8,1.0];
%     phase=[0.4,0.7,1.0];
%     for j=2:4
%         theta_r_tmp=zeros(m);
%         for i=1:m
%             theta_r_tmp(i)=atan2(wave_func((i)/m,phase(j-1))-wave_func((i-1)/m,phase(j-1)),1/m)*180/pi;
%         end
%         theta_r(j,1)=theta_r_tmp(1);
%         for i=2:m
%             theta_r(j,i)=theta_r_tmp(i)-theta_r_tmp(i-1);
%         end
%     end
%     
%     
%     
%     
% theta_r = theta_r * pi/180;%[rad]
% 
% %%
% 
% % 
% pp = zeros(2,m+2);
% pr = zeros(2,m+1);
% 
% pp(1,1) = 0;
% pp(2,1) = 0;
% 
% theta_r2 = theta_r(num,:);
% 
% theta(1,1) = theta_r2(1,1);
% 
% for i = 2:m
%     theta(1,i) = theta(1,i-1) + theta_r2(1,i);
% end
% 
%    
%     for i = 2:m+1
%         pp(1,i) = l(1,i-1) * cos(theta(1,i-1)) + pp(1,i-1);
%         pp(2,i) = l(1,i-1) * sin(theta(1,i-1)) + pp(2,i-1);
%         
%         pr(1,i) = l(1,i-1) * cos(theta(1,i-1)) + pr(1,i-1);
%         pr(2,i) = l(1,i-1) * sin(theta(1,i-1)) + pr(2,i-1);
%     end
% 
% x_1 = zeros(1,m+2);
% y_1 = zeros(1,m+2);
% 
% x_1(1,1) = -l(1,1);
% y_1(1,1) = 0;
% 
% for i=1:m+1
%     x_1(1,i+1) = pp(1,i);
%     y_1(1,i+1) = pp(2,i);
% end
% 
% 
% 
% tt = linspace(0,2*pi,100);
% 
% po1 = zeros(2,m+1);
% po2 = zeros(2,m+1);
% 
% po1_x = zeros(1,m);
% po1_y = zeros(1,m);
% 
% po2_x = zeros(1,m);
% po2_y = zeros(1,m);
% 
% %lh = 5;
% 
% for i=1:m
%     
%     po1_x(1,i)= lh * cos(theta(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%     po1_y(1,i)= lh * sin(theta(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
%     
%     po2_x(1,i)= -1 * lh * cos(theta(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%     po2_y(1,i)= -1 * lh * sin(theta(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
%     
%     po1(1,i) = po1_x(1,i);
%     po1(2,i) = po2_x(1,i);
%     
%     po2(1,i) = po1_y(1,i);
%     po2(2,i) = po2_y(1,i);
%     
% end
% 
% po1(1,m+1) = -l(1,1)/2;
% po1(2,m+1) = -l(1,1)/2;
%     
% po2(1,m+1) = lh;
% po2(2,m+1) = -1*lh;
% 
% figure
% %plot(x_1,y_1,'b-o','LineWidth',0.5,'MarkerSize',1,'MarkerFaceColor','b');
% %plot(x_1,y_1,'b-','LineWidth',2);
% % plot(x_1,y_1,'k-','LineWidth',1.5,'MarkerEdgeColor','k');
% % 
% % for i=1:m+1
% %     hold on
% %     plot(po1(:,i),po2(:,i),'k-','LineWidth',1.5);
% %     
% % end
% 
%     fig1 =plot(x_1,y_1,'-','LineWidth',0.7,'MarkerEdgeColor','k','MarkerSize',0.7);
%     fig1.Color = [184/255 184/255 184/255];
%     %plot(x_1,y_1,'r-','LineWidth',2,'MarkerEdgeColor','r','MarkerSize',2);
%     %plot(x_1,y_1,'b-o','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);
% % 
%     for i=1:m+1
%         hold on
%         fig1 = plot(po1(:,i),po2(:,i),'k-','LineWidth',0.7);
%         fig1.Color = [184/255 184/255 184/255];
%     end
% 
% %grid on
% %l_axis = l(1,1)*(m*2);
% l_axis = l(1,1)*(m*1.1);
% axis([-1*l_axis l_axis -1*l_axis l_axis]);
% pbaspect([1 1 1]);
% print(gcf,'-djpeg','-r300',fig_file2)
% %saveas(gcf,fig_file2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%