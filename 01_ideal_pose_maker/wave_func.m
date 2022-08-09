function [y] = wave_func(x,t)
%WAVE_FUNC この関数の概要をここに記述
%   詳細説明をここに記述
amp=0.3; %山の高さ
sigma=0.2; %山の幅
y=amp*exp(-((1/sigma*(x-t))^2));
end

