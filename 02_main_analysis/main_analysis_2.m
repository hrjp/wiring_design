%�听�����̓v���O����
% �ϐ����T���v�����ɑΉ�����ver
clear all

filename_1 = sprintf('output/main_matrix');
filename_2 = sprintf('output/main_Y');
filename_3 = sprintf('output/main_avg');
filename_4 = sprintf('output/main_Y2');
filename_5 = sprintf('output/main_Y0');
filename_6 = sprintf('output/main_A_hukugen');
filename_7 = sprintf('output/main_ratios_d');
filename_8 = sprintf('output/main_theta0');

%m = 100;         %�֐ߐ�
min_d = 0.7; %�ݐϊ�^��

%�ڕW�f�[�^�̃C���|�[�g
A_raw = readmatrix('../01_ideal_pose_maker/output/ideal_theta_abs.csv');
A = A_raw(2:size(A_raw,1),:);

components = size(A,2);       %�ϐ��̐�
number = size(A,1);           %�f�[�^��

%�e��v�f�̕��ςƕW���΍������߂�
avg = mean(A);
s = std(A);

%�ǂݍ��񂾃f�[�^�̊e��𕽋�0,�W���΍�1�ŕW��������
Normalization = zeros(number, components);      %�e�ʂ̊m��(�v�Z�̍������̂���)
for n = 1:number
    for comp = 1:components
        Normalization(n,comp) = ((A(n,comp) - avg(comp))/ s(comp));
    end
end

%�f�[�^�̑��֍s������߂�
R = corrcoef(A);
%���֍s��̌ŗL�l�̌v�Z
[V_a,D_a] = eig(R);     % V : �ŗL�x�N�g��, D : �ŗL�l�̍s��

%�ŗL�l�̒��o
d_a = diag(D_a);

%���ѕς� d�̑傫����
[d,ind] = sort(d_a,'descend');
D = D_a(ind,ind);
V = V_a(:,ind);

%��^��
ratios_d = (d)/sum(d);

%�݌v��^��
sum_ratios_d = zeros( size(d,1),1 );
sum_ratios_d(1) = ratios_d(1);
for i = 2:size(d,1)
    sum_ratios_d(i) = sum_ratios_d(i-1) + ratios_d(i);
end

%������
P_C_S = zeros(components,components);

% �s��̍팸
j = 0;
check = 0;
for i = 1 : components    %�s��D�̃T�C�Y��(componets * components)
    if check == 0       %���O����
        j = j + 1;
        P_C_S(:,j) = V(:,i); %Principal Components Score
        if sum_ratios_d(i) > min_d  %�ݐϊ�^����min_d�𒴂���܂Ŏg��
            check = 1;
            disp(['�ݐϊ�^�� : ',num2str(sum_ratios_d(i))]);
            disp(['����p���C���[�{�� : ',num2str(i)]);
        end
    end
end

P_C_S = P_C_S(1:components,1:j);
Y = Normalization * P_C_S;
components2 = size(Y,2);

%�m�F
N2 = P_C_S *  Y.';
che_1 = Normalization - N2.';
sum_che_1 = sum(abs(che_1),2);
%

% �W������߂�
P_C_S_2 = zeros(components,components2);
for comp2 = 1:components2
    for comp = 1:components
        P_C_S_2(comp,comp2) = s(comp)*P_C_S(comp,comp2);
    end
end

pat_n = size(P_C_S_2,2);

%�m�F
A_f = P_C_S_2 *  Y.' ;
A_fin = zeros(components,number);
for i = 1:number
    A_fin(:,i) = A_f(:,i) + avg(:);
end
che_2 = A_fin - A.';
sum_che_2 = sum(abs(che_2),1);
%

%Y�̕����I�t�Z�b�g��t����
%���������ׂĐ���Y2
[Ymin,Ymin_ind] = min(Y);
Y2 = zeros(number,components2);
Ymin_2 = zeros(number,components2);
for i = 1:number
    Y2(i,:) = Y(i,:)-Ymin(1,:);
    Ymin_2(i,:) = Ymin(1,:);
end
%�I�t�Z�b�gY0
Y0 = P_C_S_2 *  Ymin_2.';

%�m�F
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



%EXCEL�o��
% �ϊ��s�� P_C_S_2
filename = append(filename_1,'.csv');
writematrix(P_C_S_2,filename)

% �听���T���v��
filename = append(filename_2,'.csv');
WRITE = Y;
writematrix(WRITE,filename)

% A����
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

%�@����theta
filename = append(filename_6,'.csv');
WRITE = A_fin.';
writematrix(WRITE,filename)

%�@��^���@ratios_d
filename = append(filename_7,'.csv');
WRITE = ratios_d;
writematrix(WRITE,filename)

%�@theta0
filename = append(filename_8,'.csv');
WRITE = theta0;
writematrix(WRITE,filename)
