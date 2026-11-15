function [passive_tendon_length_change, new_point] = ventral_two_cable_bend_fun(D, gravity_direction)
%% logarithmic spiral parameters
b = 0.16; 
interval = pi/6; %30度
kp = 1/exp(b*interval);%縮小倍率
theta01 = 3.5 * pi;%最大節起始位置在630度
theta02 = 0.5 * pi +  interval;%最小節終點在120度
theta1 = theta01 : -interval : theta02 - 2*pi;%形成30度的等差數列，外螺旋
theta2 = theta01 : -interval : theta02;%內螺旋
k = 0.5 + 0.5 * exp (-2 * b * pi);%內外螺旋半徑比例0.68
section_number = round((max(theta2)-min(theta2))/interval + 1);%17個節點，但有18節
alpha_value = 1;

%% the radius R for the first section 
ya = exp(b * theta01) * sin(theta01);%外螺旋起點y
yb = k * exp(b * theta01) * sin(theta01);%內螺旋起點y
deltay = abs(ya - yb);
xn = -exp(b * (theta01 - interval)) * cos(theta01 - interval);
yn = exp(b * (theta01 - interval)) * sin(theta01 - interval);
degree = atan(xn / (yn - ya));
R = deltay * sin(degree);%用來計算下一段反推a


%% define contour logarithmic spiral
a = 44 / R; % make sure R = 44 mm
x1 = - a * exp (b * theta1) .* cos(theta1);
y1 = a * exp (b * theta1) .* sin(theta1);
z1 = zeros(1, length(theta1));%一個陣列，有18個0

x2 = - k * a * exp (b * theta2) .* cos(theta2);
y2 = k * a * exp (b * theta2) .* sin(theta2);
z2 = zeros(1, length(theta2));

% length of the sections
L = zeros(1,section_number);
sec1_length = sqrt(x2(2)^2 + (y2(2) - y2(1))^2 + z2(2)^2);
sitex = zeros(1,section_number);
sitex(1) = sec1_length;%這兩個相似且後面幾乎沒用到
L(1) = sec1_length;


% right side length of the first section
side1 = sqrt((x2(2) - x1(14))^2 + (y2(2) - y1(14))^2 + (z2(2) - z1(14))^2);%內螺旋第二節等於外螺旋第14節

%% positions for the threading holes
sec1_height = 44;
sec1_cable_y = sec1_height;
upper_cable_x = 19.4797;
lower_cable_x = 47.4093;
ori_coordinate = [[1, 0, 0, 0] 
                  [0, 1, 0, 0]
                  [0, 0, 1, 0] 
                  [0, 0, 0, 1]];

site1 = [upper_cable_x, sec1_cable_y, 0];
site2 = [upper_cable_x, -sec1_cable_y/2, sqrt(3) * sec1_cable_y / 2];%(上,-R/2,...)
site3 = [upper_cable_x, -sec1_cable_y/2, - sqrt(3) * sec1_cable_y / 2];

X = [lower_cable_x, sec1_cable_y, 0];
Y = [lower_cable_x, -sec1_cable_y/2, sqrt(3) * sec1_cable_y / 2];
Z = [lower_cable_x, -sec1_cable_y/2, - sqrt(3) * sec1_cable_y / 2];

T1 = ori_coordinate * trans(sec1_length, 0, 0);%trans是一個平移自訂函數
T1_point = [T1(1,4), T1(2,4), T1(3,4)];%從T1提取第四行前三個元素，T1是關節點
site4 = 6/10 * X + 4/10 * T1_point;%邊緣60%，40%關節點
site5 = 6/10 * Y + 4/10 * T1_point; 
site6 = 6/10 * Z + 4/10 * T1_point;

site1 = 6/10 * site1;%hole_point = 0.6
site2 = 6/10 * site2;
site3 = 6/10 * site3;


% beam theory parameters
%% bending moment I (stiffness)
radius = zeros(1 , (section_number - 1));
r = 3.825 * 1e-3;%彈性中心軸
for i = 1 : (section_number - 1)
    radius(i) = r * (kp ^ (i - 1));%0.68縮減比
end

stiffness = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    stiffness(i) = (pi / 4) * ( radius(i) ^ 4 );
end

%% Youngs Modulus E
E = 16.9e6;

%% moment M (Torque)
Torque = zeros(1 , (section_number - 1));
%assume a very small torque, where the smallest joint reaches its limit 30^
torque = 0.0127;
for i = 1 : (section_number - 1)
    Torque(section_number - i) = torque * kp ^ (-1/2) ^ (i - 1); %assume the torque increase as the section becomes larger
end
%kp^-1/2是經驗

%% joint length L  
L_joint = zeros(1 , (section_number - 1));
L_init = 2.10 * 1e-3;%關節起始長度
for i = 1 : (section_number - 1)
    L_joint(i) = L_init * (kp ^ (i - 1));
end

%% changing angles: theta_bend, and end_length
theta_bend = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    theta_bend(i) = Torque(i) * L_joint(i) / (E * stiffness(i));
end

% r_hole: the length from joint to threading hole 穿線洞
r_hole = zeros(1 , (section_number - 1));
r_hole_init = 28.87 * 1e-3;
for i = 1 : (section_number - 1)
    r_hole(i) = r_hole_init * kp ^ i;
end

end_length = zeros(1 , (section_number - 1));%end_length: the cable length after rotation
for i = 1 : (section_number - 1)
    end_length(i) = sqrt(2 * r_hole(i) ^ 2 * (1 - cos(pi / 6 + theta_bend(i))));
end

%% original cable length between sections (gravity 90°)
delta_length = zeros(1 , (section_number - 1));
origin_length = zeros(1 , (section_number - 1));
origin_length_init = 13.74 * 1e-3;%最大節起始線長
for i = 1 : (section_number - 1)
    origin_length(i) = origin_length_init * kp ^ (i - 1);
end

for i = 1 : (section_number - 1)
    delta_length(i) = end_length(i) - origin_length(i);
end

sum_origin_length = sum(origin_length);%全長

%% distribution portion 每一節分到的伸長量
Portion = zeros(1 , (section_number - 1));
for i = 1 : (section_number - 1)
    Portion(i) = delta_length(i) / sum(delta_length);
end

%% maximum cable length between sections: t_bend
hole_point = 0.6;
t = zeros(1, (section_number - 1));%沒什麼作用的變數
t_bend = zeros(1, (section_number - 1));
a = hole_point * side1;
b = hole_point * side1;

t(section_number - 1) = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cos(pi/6));
t_bend(1) = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cos(pi/3)) * 1e-3;

t_bend_sum = sum(t_bend);

for i = 2 : (section_number - 1) 
    t_bend(i) = t_bend(i - 1) * kp;
    t_bend_sum = t_bend_sum + t_bend(i);
end

%% initial cable lengths under certain gravity direction
gravity_angle = solve_angle(deg2rad(gravity_direction));%呼叫solve
gravity_tendon_length = zeros(1,section_number - 1);
for i = 1 : section_number - 1
    a = r_hole(i);
    b = r_hole(i);
    gravity_tendon_length(i) = angle2tendon(a, b, pi/6 + deg2rad(gravity_angle(i)));
end

gravity_tendon_length_sum = sum(gravity_tendon_length);

%% distribute actuating lengths to each sections
new_tendon_length = gravity_tendon_length;
distribute_length = zeros(1, (section_number - 1));
for i = (section_number - 1) : -1 : 1
    distribute_length(i) = distribute_length(i) + D * 1e-3 * Portion(i);
    if (gravity_tendon_length(i) + distribute_length(i)) < t_bend(i)
        new_tendon_length(i) = new_tendon_length(i) + distribute_length(i);
    else 
        new_tendon_length(i) = t_bend(i);
        addtion = distribute_length(i) - (t_bend(i) - origin_length(i));%外溢
        sum_rest = sum(delta_length(1 : i - 1));
        for j = 1 : i - 1
            distribute_length(j) = distribute_length(j) + addtion * delta_length(j) / sum_rest;%加權轉移
        end
    end
end

final_angle = zeros(1, (section_number - 1));
for k = 1 : (section_number - 1)
    a = r_hole(k);
    b = r_hole(k);
    final_angle(k) = tendon2angle(a, b, new_tendon_length(k));%tendon2angle(a, b, L)：餘弦定理的反函數
end

%% calculate joint coordinate

rotate_angle = zeros(1, (section_number - 1));
for k = 1 : (section_number - 1)
    rotate_angle(k) = 0.5236 - final_angle(k);%0.5236 弧度：即pi/6
end

new_coordinate = zeros(section_number, 4, 4);% 儲存 17 個關節各自的 4x4 齊次變換矩陣
new_point = zeros(section_number, 4);% 儲存 17 個關節中心點的 3D 齊次座標 [x, y, z, 1]

last_point1 = zeros(section_number - 1, 4);% 關節盤上第 1 個特徵參考點 3D 齊次座標
last_point2 = zeros(section_number - 1, 4);
last_point3 = zeros(section_number - 1, 4);

active_point = zeros(section_number - 1, 4);% 主動拉繩 (Active Cable) 穿線孔的 3D 空間座標
passive_point1 = zeros(section_number - 1, 4);
passive_point2 = zeros(section_number - 1, 4);

passive_tendon_length = zeros(section_number - 1,1);% 儲存最終被動繩長度
active_tendon_length = zeros(section_number - 1,1);


%% drawing
figure();
for k = 1 : (section_number)
    if k == 1 % translate and rotate with respect to the origin
        ori_coordinate = [[1, 0, 0, 0] 
                          [0, 1, 0, 0]
                          [0, 0, 1, 0] 
                          [0, 0, 0, 1]];

        new_coordinate(k,:,:) = ori_coordinate * rotate_z(deg2rad(-gravity_direction)) * trans(sec1_length, 0, 0) *  rotate_z(rotate_angle(k)) ; %轉移矩陣
        new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]';%提取座標


        site4_ = site4 * rotate_z_(deg2rad(gravity_direction));
        site5_ = site5 * rotate_z_(deg2rad(gravity_direction));
        site6_ = site6 * rotate_z_(deg2rad(gravity_direction));

        scatter3(site4_(1), site4_(2), site4_(3), 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(site5_(1), site5_(2), site5_(3), 'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(site6_(1), site6_(2), site6_(3), 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;


        last_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site4 * (kp ^ k), 1]';
        last_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site5 * (kp ^ k), 1]';
        last_point3(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site6 * (kp ^ k), 1]';

        scatter3(last_point1(k,1), last_point1(k,2), last_point1(k,3), 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point2(k,1), last_point2(k,2), last_point2(k,3), 'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point3(k,1), last_point3(k,2), last_point3(k,3), 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;


        active_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site1 * (kp ^ k), 1]';
        passive_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site2 * (kp ^ k), 1]';
        passive_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site3 * (kp ^ k), 1]';
               
        scatter3(active_point(k,1), active_point(k,2), active_point(k,3), 'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;%Alpha 透明度
        scatter3(passive_point1(k,1), passive_point1(k,2), passive_point1(k,3), 'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[1 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point2(k,1), passive_point2(k,2), passive_point2(k,3), 'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;


        plot3([site4_(1), active_point(k,1)], [site4_(2), active_point(k,2)], [site4_(3), active_point(k,3)]);hold on;
        plot3([site5_(1), passive_point1(k,1)], [site5_(2), passive_point1(k,2)], [site5_(3), passive_point1(k,3)]);hold on;
        % plot3([site6(1), passive_point2(k,1)], [site6(2), passive_point2(k,2)], [site6(3), passive_point2(k,3)]);hold on;

        passive_tendon_length(k) = dis(site5_, passive_point1(k,:));
        active_tendon_length(k) = dis(site4_, active_point(k,:));
       
    elseif k < section_number % translate and rotate with respect to the last coordinate
        last_coordinate = reshape(new_coordinate(k - 1,:,:), [4,4]);
        new_coordinate(k,:,:) = last_coordinate * trans(sec1_length * kp ^ (k - 1), 0, 0) *  rotate_z(rotate_angle(k));
        new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]'; 
       

        last_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site4 * (kp ^ k), 1]';
        last_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site5 * (kp ^ k), 1]';
        last_point3(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site6 * (kp ^ k), 1]';

        scatter3(last_point1(k,1), last_point1(k,2), last_point1(k,3), 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point2(k,1), last_point2(k,2), last_point2(k,3), 'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(last_point3(k,1), last_point3(k,2), last_point3(k,3), 'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
    
    
        active_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site1 * (kp ^ k), 1]';
        passive_point1(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site2 * (kp ^ k), 1]';
        passive_point2(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [site3 * (kp ^ k), 1]';
        

        scatter3(active_point(k,1), active_point(k,2), active_point(k,3), 'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point1(k,1), passive_point1(k,2), passive_point1(k,3), 'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[1 1 0], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;
        scatter3(passive_point2(k,1), passive_point2(k,2), passive_point2(k,3), 'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1], 'MarkerFaceAlpha', alpha_value,'MarkerEdgeAlpha', alpha_value);hold on;



        plot3([last_point1(k-1,1), active_point(k,1)], [last_point1(k-1,2), active_point(k,2)], [last_point1(k-1,3), active_point(k,3)]);hold on;
        plot3([last_point2(k-1,1), passive_point1(k,1)], [last_point2(k-1,2), passive_point1(k,2)], [last_point2(k-1,3), passive_point1(k,3)]);hold on;
        % plot3([last_point3(k-1,1), passive_point2(k,1)], [last_point3(k-1,2), passive_point2(k,2)], [last_point3(k-1,3), passive_point2(k,3)]);hold on;
        
        passive_tendon_length(k) = dis(last_point2(k - 1,:), passive_point1(k,:));
        active_tendon_length(k) = dis(last_point1(k - 1,:), active_point(k,:));
    else
            last_coordinate = reshape(new_coordinate(k - 1,:,:), [4,4]);
            new_coordinate(k,:,:) = last_coordinate * trans(sec1_length * kp ^ (k - 1), 0, 0);
            new_point(k,:) = reshape(new_coordinate(k,:,:),[4,4]) * [0,0,0,1]';

    end
end

t_rest = sum_origin_length * 1000;
t_intpu_ini = 37.6085;
t_intpu = t_intpu_ini;
for i = 2 :section_number
    t_intpu = t_intpu + t_intpu_ini * kp ^ (i - 1);
end

%% mujoco data mujoco是一個物理模擬引擎
passive_tendon_length_sum = sum(passive_tendon_length);
passive_tendon_length_change = passive_tendon_length_sum - t_rest;
mujoco_passive = t_intpu + passive_tendon_length_sum

active_tendon_length_sum = sum(active_tendon_length);
active_tendon_length_change = t_rest - active_tendon_length_sum;
mujoco_active = t_intpu + active_tendon_length_sum

new_coordinate = zeros((section_number - 1), 4, 4);
new_pointC = zeros((section_number - 1), 4);
d1 = zeros((section_number - 1), 1);

plot3(new_point(:,1)', new_point(:,2)', new_point(:,3)');
hold on;
scatter3(new_point(:,1)', new_point(:,2)', new_point(:,3)');

xlabel('x');
ylabel('y');
zlabel('z');
axis equal;


%三維平移矩陣
function T_t = trans(x, y, z)
    T_t = [[1, 0, 0, x]
           [0, 1, 0, y]
           [0, 0, 1, z]
           [0, 0, 0, 1]];
end

%繞y軸旋轉
function T_y = rotate_y(angle)
    T_y = [
        [cos(angle), 0, sin(angle)]
        [0,1,0]
        [-sin(angle), 0, cos(angle)]
        ];
end

%繞x軸旋轉
function T_x = rotate_x(angle)
    T_x = [
    [1,0 0]
    [0, cos(angle), -sin(angle)]
    [0, sin(angle), cos(angle)]
    ];
end

%繞z軸旋轉
function T_z = rotate_z(angle)
    T_z = [
    [cos(angle),  -sin(angle), 0, 0]
    [sin(angle), cos(angle), 0, 0]
    [0, 0, 1, 0]
    [0, 0, 0, 1]];
end

%已知繩長，求角度
function angle = tendon2angle(a, b, tendon_length)
    if tendon_length == 0
        angle = 0;
    else
        cos_angle = (a ^ 2 + b ^ 2 - tendon_length ^ 2) / (2 * a * b);
        angle = acos(cos_angle);
    end
end


function T_t = ptrans(x, y, z)
    T_t = [x, y, z];
end

function distance = dis(p1, p2)
    distance = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
end

%已知角度，求繩長
function tendon = angle2tendon(a, b, angle)
    tendon = sqrt(a ^ 2 + b ^ 2 - 2 * a * b * cos(angle));
end

%繞z軸旋轉
function T_z = rotate_z_(angle)
    T_z = [
    [cos(angle),  -sin(angle), 0]
    [sin(angle), cos(angle), 0]
    [0, 0, 1]
    ];
end
end