function [y] = wave_func(x,t)
%WAVE_FUNC ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
amp=0.3; %�R�̍���
sigma=0.2; %�R�̕�
y=amp*exp(-((1/sigma*(x-t))^2));
end

