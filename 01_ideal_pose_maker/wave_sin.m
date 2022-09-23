function [y] = wave_exp(x,t)
%WAVE_FUNC ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
amp=0.2; %�R�̍���
omega=6; %�p���g��
if x < (2*pi*t+pi/2)/omega
    y=amp/2*(-cos(pi*omega/(2*pi*t+pi/2)*x)+1);
else
    y=amp*sin(omega*x-2*pi*t);
end

end