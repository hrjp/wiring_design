%�����N���C���[�w���f���@���ÓI���
%�����̂�
%2�{���C���[�@���ÓI����
%�I�t�Z�b�g����

clear all

tic
%�ݒ�l
%str1 ="0125";
str1 ="0208011";

%�ݒ�l
m = 50;             %�����N��
n = m;
l = zeros(1,m);     %�����N���� [mm]
l(1,:) = 6;
%l(1,50) = 3;

%�f�[�^�̃C���|�[�g
M = readmatrix('../02_main_analysis/output/main_matrix_abs_s.csv');
Y_theta0 = readmatrix('../02_main_analysis/output/main_theta0_abs_s.csv');
Y2 = readmatrix('../02_main_analysis/output/main_Y2.csv');
%Y0 = readmatrix('main_Y0.csv');
%avg = readmatrix('main_avg.csv');
T_pat = readmatrix('output/pasig_T.csv');
pasig_pat = readmatrix('output/pasig_loop2.csv');

pat_n = size(M,2);      %�听���v�f��
%pat = 0;                %�p�^�[�����@�C���f�b�N�X

%�ڕW�l
T = zeros(pat_n,1);                      %����

T_r1 = zeros(pat_n,1);
T_s1 = zeros(pat_n,1);
dT1 = zeros(pat_n,1);

% �ǂ̃p�^�[�����Č�����̂�
%number_i = 4;          %�g�p�p�^�[��(pay_i < pat_n)

number_n = size(Y2,1);

for number_i = 4:4
    
    %�ڕW���́@���
    for pat = 1:pat_n
        %T_r1�̑��
        T_r1(pat,1) = T_pat(number_i,pat);
        T_s1(pat,1) = 0;
    end

    [T_max,T_max_i] = max(T_r1);
    % T_r1�̂����C�ő�̃��m��dT=0.01
    dT1(T_max_i,1) = 0.01;
    % ����T_r1_max����Ƃ��Č��߂�D
    for pat = 1:pat_n
        if pat == T_max_i
        else
            %dT1(pat,1) = dT1(T_max_i,1) * T_r1(pat,1)/T_r1(T_max_i,1);
            dT1(pat,1) =0;
        end
    end

    % T_r1�̂����C�ő�̃��m����{�Ƃ��Đi�߂�
    T_r = T_r1(T_max_i,1);
    T_s = T_s1(T_max_i,1);
    dT = dT1(T_max_i,1);

    lp_kari = (T_r-T_s)/dT;           %�v�Z���[�v��
    lp = round(lp_kari);

    T_r1(1,1) = T_pat(number_i,1);

    %�I�t�Z�b�g�p�x ���
    theta0 = zeros(1,size(Y_theta0,1));
    for i = 1:size(Y_theta0,1)
        theta0(1,i) = Y_theta0(i,number_i);
        %theta0(1,i) = Y0(i,number_i) + avg(1,i);
    end
    %�I�t�Z�b�g�p�x ��Όn
    theta0_abs = zeros(1,size(Y_theta0,1));
    for i = 1:size(Y_theta0,1)
        if i == 1
            theta0_abs(1,i) = theta0(1,i);
        else
            theta0_abs(1,i) = theta0_abs(1,i-1)+theta0(1,i);
        end
    end

    loop_fin = 1;
    % %loop�L�^�p
    theta_r_loop = zeros(loop_fin,m);
    theta_fin_loop = zeros(loop_fin,m);
    theta_gosa_loop = zeros(loop_fin,m);

    pasig_loop = zeros(2,m,loop_fin);
    pin_exist_loop = zeros(loop_fin,n);
    sum_theta_gosa_loop = zeros(loop_fin,1);
    WRITE = zeros(5*loop_fin,m+2);
    WRITE2 = zeros(loop_fin,2); 
    WRITE3 = zeros(loop_fin,n+1); 

    filename1 = sprintf('number_%d_',number_i);
    fig_file = append('output/pin_fig_',filename1,str1,'.jpeg');
    fig_file2 = append('output/pin_fig0_',filename1,str1,'.jpeg');
    video_file = append('output/pin_movie_',filename1,str1,'.avi');
    
    l_pa = zeros(2,m);  %pa_1-p_1�̋����x�N�g���@l_pa(1,m)=�����N���s����, l_pa(2,m)=�����N�������� : 1����m��m��

    l_pa(1,:) = l(1,1)/2;
    l_pa(2,:) = 0;
    l_pa0 = zeros(2,pat_n);
    for pat = 1:pat_n
        l_pa0(1,pat) = -1*l(1,1)/2;
        %l_pa0(1,pat) = pasig_pat(3*pat-1,1);
        l_pa0(2,pat) = pasig_pat(3*pat,1);
    end

    pasig_min = 10^(-4);

    % pasig �̑��
    pasig = zeros(2,n,pat_n);
    for i = 1:pat_n
        for j = 1:n
            %pasig(1,j,i) = l(1,j)/2;
            pasig(1,j,i) = pasig_pat(3*i-1,j+1);
            pasig(2,j,i) = pasig_pat(3*i,j+1);
        end
    end

    %�ϐ�
    p = zeros(2,m+1);     %�����N�֐߈ʒu�x�N�g�� 1=x, 2=y
    theta = zeros(1,m); %�����N�֐ߊp�x[rad]
    thetasig = zeros(1,m); %�����N�֐ߊp�x[rad]

    pa = zeros(2,n+1,pat_n);    %�����N���C���[�ʂ����S�x�N�g�� 1=x, 2=y : 0����m��m+�P��
    eb = zeros(2,n+1,pat_n);     %����tb�P�ʃx�N�g�� : 0����m��m+�P��

    q = zeros(2,m,n,pat_n);
    %q = zeros(2,m,pat_n);

    p_r = zeros(2,m+1);     %�ڕW�����N�֐߈ʒu�x�N�g�� 1=x, 2=y
    pa_r = zeros(2,n);    %�ڕW�����N���C���[�ʂ����S�x�N�g�� 1=x, 2=y : 0����m��m+�P��
    eb_r = zeros(2,n);     %�ڕW����tb�P�ʃx�N�g�� : 0����m��m+�P��
    ebsig_r = zeros(2,n);     %�ڕW����tb�P�ʃx�N�g�� : 0����m��m+�P��

    %k�v�Z�p
    Q = zeros(m,m);         %�v�Z�p�@�֐ߊp�x����s��
    ABe = zeros(m,pat_n);       %�v�Z�p�@���͍�����s��
    u = zeros(m,1);       %�v�Z�p�@T������
    u1 = zeros(m,pat_n);
    L = zeros(1,pat_n);

    %���ÓI�ߒ��v�Z�p
    dtheta = zeros(1,m); %�����N�֐ߊp�x�ω�[rad]
    K = zeros(m,m);     %�΂˃g���N�W����
    K0 = zeros(m,m);     %�I�t�Z�b�g���΂˃g���N�W����
    k = zeros(1,m);
    k(1,:) = 830;

    %�L�^�l
    p_result = zeros(lp+1,2,m+1);     %�����N�֐߈ʒu�x�N�g��
    theta_result = zeros(lp+1,m); %�����N�֐ߊp�x[rad]
    pa_result = zeros(lp+1,2,n+1,pat_n);    %�����N���C���[�ʂ����S�x�N�g��
    dtheta_result = zeros(lp+1,m);
    g_kakuninn2 = zeros(lp+1,m);
    theta_r_result = zeros(lp+1,m);


    %%plot�Q���z�ŏI���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %���ÓI�ߒ��@�v�Z��lp 
    %�ȉ��C��΍��W�n
    %theta = zeros(1,m); 

    %K(1), K(m)�̂ݕʌv�Z
    K(1,1) = k(1,1) + k(1,2);
    K(1,2) = -1*k(1,2);
    K(m,m-1) = -1*k(1,m);
    K(m,m) = k(1,m);

    for i = 2:m-1
        K(i,i-1) = -1*k(1,i);
        K(i,i) = k(1,i) + k(1,i+1);
        K(i,i+1) = -1*k(1,i+1);
    end
    
    %�I�t�Z�b�g�p�x���̂΂˒萔
    K0(m,m) = k(1,m);
    for i = 1:m-1
        K0(i,i) = k(1,i);
        K0(i,i+1) = -1*k(1,i+1);
    end
    
    for i = 1:m

        theta(1,i) = theta0_abs(1,i);

    end

    %���ÓI�ߒ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %T������
    for pat = 1:pat_n
        T(pat,1) = T_s1(pat,1);
    end

    %p = zeros(2,m+1);
    for time = 1:lp+1
        % p, pa �X�V�i���ÓI�j%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %��ꃊ���N�֐߂͌Œ�
        p(1,1) = 0;
        p(2,1) = 0;
        %��񃊃��N�ȍ~�̍X�V
        for i = 2:m+1
            p(1,i) = l(1,i-1) * cos(theta(1,i-1)) + p(1,i-1);
            p(2,i) = l(1,i-1) * sin(theta(1,i-1)) + p(2,i-1);
        end

        for pat = 1:pat_n
            %pa �X�V
            pa(1,1,pat) = p(1,1) + l_pa0(1,pat);
            pa(2,1,pat) = p(2,1) + l_pa0(2,pat);
            for j = 1:n
                pa(1,j+1,pat) = p(1,j) + pasig(1,j,pat)*cos(theta(1,j)) - pasig(2,j,pat)*sin(theta(1,j));
                pa(2,j+1,pat) = p(2,j) + pasig(1,j,pat)*sin(theta(1,j)) + pasig(2,j,pat)*cos(theta(1,j));
            end

            %eb�@�X�V
            for j = 1:n
                eb(1,j,pat) = (pa(1,j,pat) -pa(1,j+1,pat)) / sqrt( (pa(1,j,pat) -pa(1,j+1,pat))^2 + (pa(2,j,pat) -pa(2,j+1,pat))^2 );
                eb(2,j,pat) = (pa(2,j,pat) -pa(2,j+1,pat)) / sqrt( (pa(1,j,pat) -pa(1,j+1,pat))^2 + (pa(2,j,pat) -pa(2,j+1,pat))^2 );
            end

        end

        %g�̌v�Z%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %q(2,m)
        q = zeros(2,m,n,pat_n);
        %q = zeros(2,m,pat_n);
        u = zeros(m,1);
        u1 = zeros(m,pat_n);

        for pat = 1:pat_n

            for i = 1:m-1
                q(1,i,i,pat) = pa(1,i+1,pat) - p(1,i);
                q(2,i,i,pat) = pa(2,i+1,pat) - p(2,i);

                %�ύX�ӏ��@pa
                q(1,i,i+1,pat) = pa(1,i+2,pat) - p(1,i+1);
                q(2,i,i+1,pat) = pa(2,i+2,pat) - p(2,i+1);

                %q(1,i,pat) = pa(1,i+1,pat) - p(1,i);
                %q(2,i,pat) = pa(2,i+1,pat) - p(2,i);
            end
            %i=m,j=n�͗�O            
            q(1,m,n,pat) = pa(1,n+1,pat) - p(1,m);
            q(2,m,n,pat) = pa(2,n+1,pat) - p(2,m);

            %q(1,m,pat) = pa(1,n+1,pat) - p(1,m);
            %q(2,m,pat) = (2,n+1,pat) - p(2,m);

            %�s��ABe ���͉^����
            for i = 1:m-1
                ABe(i,pat) = ( q(1,i,i,pat)*eb(2,i,pat) - q(2,i,i,pat)*eb(1,i,pat) ) - ( q(1,i,i+1,pat)*eb(2,i+1,pat) - q(2,i,i+1,pat)*eb(1,i+1,pat) );
                %ABe(i,pat) = ( q(1,i,pat)*eb(2,i,pat) - q(2,i,pat)*eb(1,i,pat) ) - ( q(1,i+1,pat)*eb(2,i+1,pat) - q(2,i+1,pat)*eb(1,i+1,pat) );

            end

            %ABe(1,m)�͕ʌv�Z
            ABe(m,pat) = q(1,m,n,pat)*eb(2,n,pat) - q(2,m,n,pat)*eb(1,n,pat);
            %ABe(m,pat) = q(1,m,pat)*eb(2,n,pat) - q(2,m,pat)*eb(1,n,pat);

            u1(:,pat) = ABe(:,pat);

            %T������
            u = u + u1(:,pat);
        end

        %D�̋���
        D = zeros(m,m);
        %�����v�Z(��͓I)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %���ӁF��̏��ÓI�v�Z�ɑg�ݍ���ł��悢
        %dp(i)/dtheta(i)
        d_p_theta_x = zeros(m,m);
        d_p_theta_y = zeros(m,m);

        %dpa(i+1)/dtheta(i)
        d_pa_theta_x = zeros(n,m,pat_n);
        d_pa_theta_y = zeros(n,m,pat_n);

        %d_L_theta = zeros(m,m);
        d_L_theta = zeros(n,m,pat_n);

        %deb(i+1)/dtheta(i)
        d_eb_theta_x = zeros(n,m,pat_n);
        d_eb_theta_y = zeros(n,m,pat_n);

        d_q_theta_x = zeros(m,n,m,pat_n);
        d_q_theta_y = zeros(m,n,m,pat_n);

        %d_q_theta_x = zeros(m,m,pat_n);
        %d_q_theta_y = zeros(m,m,pat_n);

        % i : p()�̃i���o�����O m
        % j : pa()�̃i���o�����O n

        for i = 1:m
            for i_theta = 1:m
                %dp
                if(i == 1)
                    d_p_theta_x(i,i_theta) = 0;
                    d_p_theta_y(i,i_theta) = 0;
                else
                    if ((i-1) < i_theta)
                        d_p_theta_x(i,i_theta) = 0;
                        d_p_theta_y(i,i_theta) = 0;
                    elseif ((i-1) > i_theta)
                        d_p_theta_x(i,i_theta) = -1 * l(1,i_theta) * sin(theta(1,i_theta));
                        d_p_theta_y(i,i_theta) = l(1,i_theta) * cos(theta(1,i_theta));
                    else
                        d_p_theta_x(i,i_theta) = -1 * l(1,i-1) * sin(theta(1,i-1));
                        d_p_theta_y(i,i_theta) = l(1,i-1) * cos(theta(1,i-1));
                    end
                end

                for pat = 1:pat_n
                    %dpa
                    if (i_theta==i)
                        d_pa_theta_x(i,i_theta,pat) = d_p_theta_x(i,i_theta) - pasig(1,i,pat)*sin(theta(1,i)) - pasig(2,i,pat)*cos(theta(1,i));
                        d_pa_theta_y(i,i_theta,pat) = d_p_theta_y(i,i_theta) + pasig(1,i,pat)*cos(theta(1,i)) - pasig(2,i,pat)*sin(theta(1,i));
                    else
                        d_pa_theta_x(i,i_theta,pat) = d_p_theta_x(i,i_theta);
                        d_pa_theta_y(i,i_theta,pat) = d_p_theta_y(i,i_theta);
                    end

                    L(1,pat) = sqrt( ( pa(1,i,pat)-pa(1,i+1,pat) )^2 + ( pa(2,i,pat)-pa(2,i+1,pat) )^2 );

                    %dL/dtheta
                    %deb/dtheta
                    if(i == 1)
                        d_L_theta(i,i_theta,pat) = ( ((pa(1,i,pat)-pa(1,i+1,pat))^2 + (pa(2,i,pat)-pa(2,i+1,pat))^2 )^(-1/2) ) * ( (pa(1,i+1,pat)-pa(1,i,pat))*(d_pa_theta_x(i,i_theta,pat) - 0) + (pa(2,i+1,pat)-pa(2,i,pat))*(d_pa_theta_y(i,i_theta,pat) - 0) );

                        d_eb_theta_x(i,i_theta,pat) = -1/(L(1,pat)^2) * d_L_theta(i,i_theta,pat) * (pa(1,i,pat)-pa(1,i+1,pat)) + 1/L(1,pat) * ( 0 - d_pa_theta_x(i,i_theta,pat) );
                        d_eb_theta_y(i,i_theta,pat) = -1/(L(1,pat)^2) * d_L_theta(i,i_theta,pat) * (pa(2,i,pat)-pa(2,i+1,pat)) + 1/L(1,pat) * ( 0 - d_pa_theta_y(i,i_theta,pat) );

                    else
                        d_L_theta(i,i_theta,pat) = ( ((pa(1,i,pat)-pa(1,i+1,pat))^2 + (pa(2,i,pat)-pa(2,i+1,pat))^2 )^(-1/2) ) * ( (pa(1,i+1,pat)-pa(1,i,pat))*(d_pa_theta_x(i,i_theta,pat) - d_pa_theta_x(i-1,i_theta,pat)) + (pa(2,i+1,pat)-pa(2,i,pat))*(d_pa_theta_y(i,i_theta,pat) - d_pa_theta_y(i-1,i_theta,pat)) );

                        d_eb_theta_x(i,i_theta,pat) = -1/(L(1,pat)^2) * d_L_theta(i,i_theta,pat) * (pa(1,i,pat)-pa(1,i+1,pat)) + 1/L(1,pat) * ( d_pa_theta_x(i-1,i_theta,pat) - d_pa_theta_x(i,i_theta,pat) );
                        d_eb_theta_y(i,i_theta,pat) = -1/(L(1,pat)^2) * d_L_theta(i,i_theta,pat) * (pa(2,i,pat)-pa(2,i+1,pat)) + 1/L(1,pat) * ( d_pa_theta_y(i-1,i_theta,pat) - d_pa_theta_y(i,i_theta,pat) );

                    end
                end    

            end

        end

        for pat = 1:pat_n
            %dq_xy(m,n,m)
            for i = 1:m-1
                for i_theta = 1:m
                    d_q_theta_x(i,i,i_theta,pat) = d_pa_theta_x(i,i_theta,pat) - d_p_theta_x(i,i_theta);
                    d_q_theta_y(i,i,i_theta,pat) = d_pa_theta_y(i,i_theta,pat) - d_p_theta_y(i,i_theta);

                    d_q_theta_x(i,i+1,i_theta,pat) = d_pa_theta_x(i+1,i_theta,pat) - d_p_theta_x(i+1,i_theta);
                    d_q_theta_y(i,i+1,i_theta,pat) = d_pa_theta_y(i+1,i_theta,pat) - d_p_theta_y(i+1,i_theta);

                    %d_q_theta_x(i,i_theta,pat) = d_pa_theta_x(i,i_theta,pat) - d_p_theta_x(i,i_theta,pat);
                    %d_q_theta_y(i,i_theta,pat) = d_pa_theta_y(i,i_theta,pat) - d_p_theta_y(i,i_theta,pat);

                end
            end
            %i=m,j=n�͗�O
            for i_theta = 1:m
                d_q_theta_x(m,n,i_theta,pat) = d_pa_theta_x(n,i_theta,pat) - d_p_theta_x(m,i_theta);
                d_q_theta_y(m,n,i_theta,pat) = d_pa_theta_y(n,i_theta,pat) - d_p_theta_y(m,i_theta);       

                %d_q_theta_x(m,i_theta,pat) = d_pa_theta_x(n,i_theta,pat) - d_p_theta_x(m,i_theta,pat);
                %d_q_theta_y(m,i_theta,pat) = d_pa_theta_y(n,i_theta,pat) - d_p_theta_y(m,i_theta,pat);       

            end
        end

        %�����v�Z�I��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        D_1 = zeros(m,m,pat_n);

        for i = 1:m
            for i_theta = 1:m
                for pat = 1:pat_n
                    if (i ==m)
                        D_1(i,i_theta,pat) = ( ( d_q_theta_x(i,i,i_theta,pat)*eb(2,i,pat) + q(1,i,i,pat)*d_eb_theta_y(i,i_theta,pat) ) - ( d_q_theta_y(i,i,i_theta,pat)*eb(1,i,pat) + q(2,i,i,pat)*d_eb_theta_x(i,i_theta,pat) ) );
                        %D_1(i,i_theta,pat) = ( ( d_q_theta_x(i,i_theta,pat)*eb(2,i,pat) + q(1,i,pat)*d_eb_theta_y(i,i_theta,pat) ) - ( d_q_theta_y(i,i_theta,pat)*eb(1,i,pat) + q(2,i,pat)*d_eb_theta_x(i,i_theta,pat) ) );

                    else
                        D_1(i,i_theta,pat) = ( ( d_q_theta_x(i,i,i_theta,pat)*eb(2,i,pat) + q(1,i,i,pat)*d_eb_theta_y(i,i_theta,pat) ) - ( d_q_theta_y(i,i,i_theta,pat)*eb(1,i,pat) + q(2,i,i,pat)*d_eb_theta_x(i,i_theta,pat) ) ) - ( ( d_q_theta_x(i,i+1,i_theta,pat)*eb(2,i+1,pat) + q(1,i,i+1,pat)*d_eb_theta_y(i+1,i_theta,pat) ) - ( d_q_theta_y(i,i+1,i_theta,pat)*eb(1,i+1,pat) + q(2,i,i+1,pat)*d_eb_theta_x(i+1,i_theta,pat) ) );
                        %D_1(i,i_theta,pat) = ( ( d_q_theta_x(i,i_theta,pat)*eb(2,i,pat) + q(1,i,pat)*d_eb_theta_y(i,i_theta,pat) ) - ( d_q_theta_y(i,i_theta,pat)*eb(1,i,pat) + q(2,i,pat)*d_eb_theta_x(i,i_theta,pat) ) ) - ( ( d_q_theta_x(i+1,i_theta,pat)*eb(2,i+1,pat) + q(1,i+1,pat)*d_eb_theta_y(i+1,i_theta,pat) ) - ( d_q_theta_y(i+1,i_theta,pat)*eb(1,i+1,pat) + q(2,i+1,pat)*d_eb_theta_x(i+1,i_theta,pat) ) );

                    end
                    D(i,i_theta) = D(i,i_theta) + T(pat,1)*D_1(i,i_theta,pat);
                end
                D(i,i_theta) = D(i,i_theta) - K(i,i_theta);

            end
        end

        %dtheta = -1* inv(D) * u * dT; 
        dtheta = 0;
        for pat = 1:pat_n
            dtheta = dtheta - dT1(pat,1) * ( D \ u1(:,pat) );
        end

        dtheta_result(time,:) = dtheta.';
        theta_result(time,:) = theta;
        p_result(time,1,:) = p(1,:);
        p_result(time,2,:) = p(2,:);

        for i = 1:m
            %�������@g_kakuninn2
            g_kakuninn2(time,i) = 0;
        end

        for pat = 1:pat_n
            TABe_kari = T(pat,1)*ABe(:,pat);
            g_kakuninn2(time,:) = g_kakuninn2(time,:) + TABe_kari.';

            pa_result(time,1,:,pat) = pa(1,:,pat);
            pa_result(time,2,:,pat) = pa(2,:,pat);
        end

        theta_r_result(time,1) = theta_result(time,1);
        for i=2:m
            theta_r_result(time,i) = theta_result(time,i) - theta_result(time,i-1);
        end
        Ktheta_kari = K*theta.' - K0*theta0.';
        g_kakuninn2(time,:) = g_kakuninn2(time,:) - Ktheta_kari.';

        %�X�V
        if(time == lp+1)
        else
            theta = theta + dtheta.';

            %�@T�@���Z
            for pat = 1:pat_n
                T(pat,1) = T(pat,1) + dT1(pat,1);
            end
        end

    end

    loop = 1;
    %loop �L�^
    theta_r_loop(loop,:) = 0;
    theta_fin_loop(loop,1) = theta(1,1);
    for i = 2:m
        theta_fin_loop(loop,i) = theta(1,i) - theta(1,i-1);
    end
    for i = 1:m
        theta_gosa_loop(loop,i) = theta_r_loop(loop,i) - theta_fin_loop(loop,i);
    end

    for i = 1:m
        sum_theta_gosa_loop(loop,1) = sum_theta_gosa_loop(loop,1) + abs( theta_gosa_loop(loop,i) );
    end

    %pasig_loop(:,:,loop) = pasig(:,:);
    %pin_exist_loop = zeros(loop_fin,m);
    for ip = 1:n
        pasig_loop(:,ip,loop) = pasig(:,ip,1);
        pin_exist_loop(loop,ip) = ip;
    end

    %�������ݗp
    WRITE( (5*loop-4) , 1 ) = loop;
    WRITE( (5*loop-4) , 2:m+1 ) = theta_r_loop(loop,:);
    WRITE( (5*loop-3) , 2:m+1 ) = theta_fin_loop(loop,:);
    WRITE( (5*loop-2) , 2:m+1 ) = theta_gosa_loop(loop,:);
    WRITE( (5*loop-2) , m+2 ) = sum_theta_gosa_loop(loop,1);
    WRITE( (5*loop-1):5*loop , 2:m+1 ) = pasig_loop(:,:,loop);

    WRITE2(loop,1) = loop;
    WRITE2(loop,2) = sum_theta_gosa_loop(loop,1);

    WRITE3(loop,1) = loop;
    for ip = 1:n
        WRITE3(loop,ip+1) = pin_exist_loop(loop,ip);
    end

    %EXCEL�o��
    %data
    %filename_1 = sprintf('loop_theta_r_fin_gosa_pasig_m%d_n%d_',m,n);
    filename_1 = sprintf('output/loopdata_');
    filename = append(filename_1,filename1,str1,'.xlsx');
    writematrix(WRITE,filename)

    %sum_�덷
    filename2_1 = sprintf('output/loop_gosa_');
    filename2 = append(filename2_1,filename1,str1,'.xlsx');
    writematrix(WRITE2,filename2)

    %pin�ʒudata
    filename3_1 = sprintf('output/loop_pin_');
    filename3 = append(filename3_1,filename1,str1,'.xlsx');
    writematrix(WRITE3,filename3)

    %theta_r_result
    filename31_1 = sprintf('output/theta_result_');
    filename31 = append(filename31_1,filename1,str1,'.xlsx');
    writematrix(theta_r_result,filename31)

    toc

    %plot�Q�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     pp = zeros(2,m+2);     %�����N�֐߈ʒu�x�N�g�� 1=x, 2=y
%     pr = zeros(2,m+1);     %�����N�֐߈ʒu�x�N�g�� 1=x, 2=y
% 
%     pp(1,1) = 0;
%     pp(2,1) = 0;
% 
%     %��񃊃��N�ȍ~�̍X�V
%     for i = 2:m+1
%         pp(1,i) = l(1,i-1) * cos(theta0_abs(1,i-1)) + pp(1,i-1);
%         pp(2,i) = l(1,i-1) * sin(theta0_abs(1,i-1)) + pp(2,i-1);
% 
%         pr(1,i) = l(1,i-1) * cos(theta0_abs(1,i-1)) + pr(1,i-1);
%         pr(2,i) = l(1,i-1) * sin(theta0_abs(1,i-1)) + pr(2,i-1);
%     end
% 
%     x_1 = zeros(1,m+2);
%     y_1 = zeros(1,m+2);
% 
%     x_1(1,1) = -1*l(1,1);
%     y_1(1,1) = 0;
% 
%     for i=1:m+1
%         x_1(1,i+1) = pp(1,i);
%         y_1(1,i+1) = pp(2,i);
%     end
% 
%     %tt = linspace(0,2*pi,100);
% 
%     po1 = zeros(2,m+1);
%     po2 = zeros(2,m+1);
% 
%     po1_x = zeros(1,m);
%     po1_y = zeros(1,m);
% 
%     po2_x = zeros(1,m);
%     po2_y = zeros(1,m);
% 
%     lh = 5;
% 
%     for i=1:m
% 
%         po1_x(1,i)= lh * cos(theta0_abs(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%         po1_y(1,i)= lh * sin(theta0_abs(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
% 
%         po2_x(1,i)= -1 * lh * cos(theta0_abs(1,i)+pi()/2) + (pp(1,i) + pp(1,i+1))/2;
%         po2_y(1,i)= -1 * lh * sin(theta0_abs(1,i)+pi()/2) + (pp(2,i) + pp(2,i+1))/2;
% 
%         po1(1,i) = po1_x(1,i);
%         po1(2,i) = po2_x(1,i);
% 
%         po2(1,i) = po1_y(1,i);
%         po2(2,i) = po2_y(1,i);
% 
%     end
% 
%     po1(1,m+1) = -l(1,1)/2;
%     po1(2,m+1) = -l(1,1)/2;
% 
%     po2(1,m+1) = lh;
%     po2(2,m+1) = -1*lh;
% 
%     figure
%     plot(x_1,y_1,'b-o','MarkerSize',4);
% 
%     for i=1:m+1
%         hold on
%         plot(po1(:,i),po2(:,i),'k-','LineWidth',2);
% 
%     end
% 
%     grid on
%     %l_axis = l(1,1)*(m*2);
%     l_axis = l(1,1)*(m+1);
%     axis([-1*l_axis l_axis -1*l_axis l_axis]);
%     pbaspect([1 1 1]);
% 
%     saveas(gcf,fig_file2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %plot�Q�ŏI���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_1 = zeros(1,m+2);
    y_1 = zeros(1,m+2);

    x_1(1,1) = -1*l(1,1);
    y_1(1,1) = 0;

    for i=1:m+1
        x_1(1,i+1) = p(1,i);
        y_1(1,i+1) = p(2,i);
    end

    x_2 = zeros(pat_n,n+1);
    y_2 = zeros(pat_n,n+1);

    for pat = 1:pat_n
        for j = 1:n+1
            x_2(pat,j) = pa(1,j,pat);
            y_2(pat,j) = pa(2,j,pat);
        end
    end

    po1 = zeros(2,m+1);
    po2 = zeros(2,m+1);

    po1_x = zeros(1,m);
    po1_y = zeros(1,m);

    po2_x = zeros(1,m);
    po2_y = zeros(1,m);

    %lh = 4;
    lh = 10;

    for i=1:m

        po1_x(1,i)= lh * cos(theta(1,i)+pi()/2) + (p(1,i) + p(1,i+1))/2;
        po1_y(1,i)= lh * sin(theta(1,i)+pi()/2) + (p(2,i) + p(2,i+1))/2;

        po2_x(1,i)= -1 * lh * cos(theta(1,i)+pi()/2) + (p(1,i) + p(1,i+1))/2;
        po2_y(1,i)= -1 * lh * sin(theta(1,i)+pi()/2) + (p(2,i) + p(2,i+1))/2;

        po1(1,i) = po1_x(1,i);
        po1(2,i) = po2_x(1,i);

        po2(1,i) = po1_y(1,i);
        po2(2,i) = po2_y(1,i);

    end
    po1(1,m+1) = -l(1,1)/2;
    po1(2,m+1) = -l(1,1)/2;

    po2(1,m+1) = lh;
    po2(2,m+1) = -1*lh;


    tt = linspace(0,2*pi,100);

    figure
    plot(x_1,y_1,'r-o','LineWidth',0.5,'MarkerSize',2,'MarkerFaceColor','r');
    %plot(x_1,y_1,'k--','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',2);
%     fig1 =plot(x_1,y_1,'-','LineWidth',0.7,'MarkerEdgeColor','k','MarkerSize',0.7);
%     fig1.Color = [184/255 184/255 184/255];
    %plot(x_1,y_1,'r-','LineWidth',2,'MarkerEdgeColor','r','MarkerSize',2);
    %plot(x_1,y_1,'b-o','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);
% 
%     for i=1:m+1
%         hold on
%         fig1 = plot(po1(:,i),po2(:,i),'k-','LineWidth',0.7);
%         fig1.Color = [184/255 184/255 184/255];
%     end
%     
% 
%     for pat = 1:1
%         hold on
%         plot(x_2(pat,:),y_2(pat,:),'r-','LineWidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',0.5);
%     end
%     for pat = 2:2
%         hold on
%         plot(x_2(pat,:),y_2(pat,:),'g-','LineWidth',0.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',0.5);
%     end

    
    %hold on
    %plot(po2_x,po2_y,'b-o');
    %hold off

    grid on
    %l_axis = l(1,1)*(m*2);
    l_axis = l(1,1)*(m*1.1);
    axis([-1*l_axis l_axis -1*l_axis l_axis]);
    pbaspect([1 1 1]);
    print(gcf,'-djpeg','-r300',fig_file)
    %saveas(gcf,fig_file)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %�ȉ�����p
    x_pp = zeros(1,m+2);
    y_pp = zeros(1,m+2);

    x_ppa = zeros(pat,n+1);
    y_ppa = zeros(pat,n+1);

    frame_num = lp+1;

    figure
    %movie
    MovieData = struct('cdata',[],'colormap',[]);
    %video_file = 'move_1215_2.avi';
    outputVideo = VideoWriter(video_file);
    outputVideo.FrameRate = lp/3;
    outputVideo.Quality = 80;
    open(outputVideo);

    for time = 1:frame_num

        x_pp(1,1) = -l(1,1);
        y_pp(1,1) = 0;

        for link = 1:m+1
            x_pp(1,link+1) = p_result(time,1,link);
            y_pp(1,link+1) = p_result(time,2,link);
        end

        for pat = 1:pat_n
            for link = 1:n+1
                x_ppa(pat,link) = pa_result(time,1,link,pat);
                y_ppa(pat,link) = pa_result(time,2,link,pat);
            end   
        end

        po1 = zeros(2,m+1);
        po2 = zeros(2,m+1);

        po1_x = zeros(1,m);
        po1_y = zeros(1,m);

        po2_x = zeros(1,m);
        po2_y = zeros(1,m);

        lh = 10;

        for i=1:m

            po1_x(1,i)= lh * cos(theta_result(time,i)+pi()/2) + (p_result(time,1,i) + p_result(time,1,i+1))/2;
            po1_y(1,i)= lh * sin(theta_result(time,i)+pi()/2) + (p_result(time,2,i) + p_result(time,2,i+1))/2;

            po2_x(1,i)= -1 * lh * cos(theta_result(time,i)+pi()/2) + (p_result(time,1,i) + p_result(time,1,i+1))/2;
            po2_y(1,i)= -1 * lh * sin(theta_result(time,i)+pi()/2) + (p_result(time,2,i) + p_result(time,2,i+1))/2;

            po1(1,i) = po1_x(1,i);
            po1(2,i) = po2_x(1,i);

            po2(1,i) = po1_y(1,i);
            po2(2,i) = po2_y(1,i);

        end
        po1(1,m+1) = -l(1,1)/2;
        po1(2,m+1) = -l(1,1)/2;

        po2(1,m+1) = lh;
        po2(2,m+1) = -1*lh;

%         if(mod(time,332)==1)
%         
%             
%             %figure
%             plot(x_pp,y_pp,'g-','LineWidth',0.5,'MarkerSize',2,'MarkerFaceColor','g');
%         
%             grid on
%             %l_axis = l(1,1)*(m*2);
%             l_axis = l(1,1)*(m*1.1);
%             axis([-1*l_axis l_axis -1*l_axis l_axis]);
%             pbaspect([1 1 1]);
%             hold on
%             %print(gcf,'-djpeg','-r300',fig_file)
%             iaiaia = 0;
%             
%         end

        plot(x_pp,y_pp,'g-','LineWidth',0.5,'MarkerSize',2,'MarkerFaceColor','r');
        %plot(x_pp,y_pp,'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',2);
        %plot(x_pp,y_pp,'k-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerSize',1.5);
% 
%         for pat = 1:pat_n
%             hold on
%             plot(x_ppa(pat,:),y_ppa(pat,:),'r-o','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
%         end

%         for i=1:m+1
%             hold on
%             plot(po1(:,i),po2(:,i),'k-','LineWidth',1.5);
% 
%         end
        

        %     fig1 =plot(x_1,y_1,'-','LineWidth',0.7,'MarkerEdgeColor','k','MarkerSize',0.7);
        %     fig1.Color = [184/255 184/255 184/255];
        %     for i=1:m+1
%     for i=1:m+1
%         hold on
%         fig1 = plot(po1(:,i),po2(:,i),'k-','LineWidth',0.7);
%         fig1.Color = [184/255 184/255 184/255];
%     end

        set(gcf,'renderer','painters')
        l_axis = l(1,1)*(m*1.1);
        axis([-1*l_axis l_axis -1*l_axis l_axis]);
        pbaspect([1 1 1]);

        grid on
        drawnow;

        MovieData (time) = getframe(gcf);
        writeVideo(outputVideo,MovieData (time));

        delete(gca)

   end
    
end