function [A, B, D, E] = terminal_calculate(theta, l, w, a, b)
% 输入第i个点的theta(1)和第i+1个点的theta(2)后, 计算该2个点所确定的矩形的4个端点坐标
% 计算对应极径
r1 = a + b * theta(1);
r2 = a + b * theta(2);

% 计算对应平面坐标
x1 = r1 * cos(theta(1));
y1 = r1 * sin(theta(1));
x2 = r2 * cos(theta(2));
y2 = r2 * sin(theta(2));

% 中心点坐标C
C = zeros(1, 2);
C(1) = 1/2 * (x1 + x2);
C(2) = 1/2 * (y1 + y2);

% 计算v向量
v = zeros(1, 2);
norm = sqrt((x1 - x2)^2 + (y1 - y2)^2); % 距离的模

v(1) = (x2 - x1) / norm; % Vx
v(2) = (y2 - y1) / norm; % Vy

% 计算法向量
n = zeros(1, 2);
n(1) = v(2); % nx
n(2) = -v(1); % ny

% 端点A
A = zeros(1, 2); % 注意这里2个分量是x坐标和y坐标
A = C - 1/2 * l * v - 1/2 * w * n;

% 端点B
B = zeros(1, 2);
B = C - 1/2 * l * v + 1/2 * w * n;

% 端点D
D = zeros(1, 2);
D = C + 1/2 * l * v + 1/2 * w * n;

% 端点E
E = zeros(1, 2);
E = C + 1/2 * l * v - 1/2 * w * n;
end