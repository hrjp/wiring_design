function [y] = wave_sin2(x,t)
%WAVE_FUNC ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
amp=0.3; %�R�̍���
omega=4; %�p���g��
if 0 < omega*(x-t) && omega*(x-t) < pi
    y=amp*sin(omega*(x-t))^2;
else
    y=0.01*x;
end

end