function [y] = wave_sin2(x,t)
%WAVE_FUNC この関数の概要をここに記述
%   詳細説明をここに記述
amp=0.3; %山の高さ
omega=4; %角周波数
if 0 < omega*(x-t) && omega*(x-t) < pi
    y=amp*sin(omega*(x-t))^2;
else
    y=0.01*x;
end

end