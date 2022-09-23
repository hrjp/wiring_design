function [y] = wave_exp(x,t)
%WAVE_FUNC この関数の概要をここに記述
%   詳細説明をここに記述
amp=0.2; %山の高さ
omega=6; %角周波数
if x < (2*pi*t+pi/2)/omega
    y=amp/2*(-cos(pi*omega/(2*pi*t+pi/2)*x)+1);
else
    y=amp*sin(omega*x-2*pi*t);
end

end