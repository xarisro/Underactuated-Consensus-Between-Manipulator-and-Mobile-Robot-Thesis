%% Online Control
%Online.m

% General Simulation and robot parameters
robot = importrobot('iiwa14.urdf');
robot.DataFormat = 'column';

end_eff_name = 'iiwa_link_ee_kuka';

inc= 0.7;
T1 = [cos(inc) -sin(inc) 0 0; sin(inc) cos(inc) 0 0; 0 0 1 0; 0 0 0 1];
T2 = [cos(inc) 0 sin(inc) 0; 0 1 0 0; -sin(inc) 0 cos(inc) 0.1; 0 0 0 1];
Ru = T1 * T2; %Incline transform
% Ru = eye(4); %No incline
setFixedTransform(robot.Bodies{1}.Joint,Ru)

% dofs = [1 2 6 7];
dofs = [2 3 4 6];
ndofs = length(dofs);
total_dofs = 7;

dt = 0.001;
Time = 20;
samples = Time / dt;

%% Manipulator Parameter initialization
q = zeros(total_dofs,samples);
q(:,1) = [ 0; 0; pi/4; -pi/6; 0; pi/8; pi/3]; %Initial configuration
qdot = zeros(total_dofs,1);

minlim = zeros(ndofs,1);
maxlim = zeros(ndofs,1); %Joint Limits

for i=1:ndofs
   minlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(1);
   maxlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(2);
end

J = zeros(6,ndofs);
Re = zeros(3,3);

%% Error initialization
ez= zeros(1,samples);
ephi = zeros(1,samples);
ep = zeros(2,samples);
emtheta = zeros(1,samples);
Hplot = zeros(1,samples);

%% Direction of apporach
trans = getTransform(robot,q(:,1),end_eff_name);
zin = trans(3,4);
zf = 0.5;
zmax = max(zf,zin) + 0.1;

T = 8;
S = T / dt;
t = linspace(0, T, S);

w = 2 * pi / T;
zd1 = 1/2 * (1 - cos(w * t)) * (zmax - zin) + zin;
zd2 = 1/2 * (1 + cos(w * t)) * (zf - zmax) + zmax;

%Desired z trajectory
zd = zd1 .* (1 - heaviside(t - T/2)) + zd2 .* heaviside(t - T/2);
zd = [zd zd(S) * ones(1, samples - S)];

%Its derivative
zddot = (1/2*w*sin(w*t)*(zmax - zin)) .* (1 - heaviside(t - T/2)) + ...
    (-1/2 * w * sin(w*t)*(zf - zmax)) .* heaviside(t - T/2);
zddot = [zddot zeros(1,samples-S)];

%% Unicycle Mobile Robot Parameter initialization
xm = zeros(3,samples);
xm(:,1) = [2.3; -2.3; 0.6]; %x,y,theta

b = 0.005; %Final Distance from the manipulator

ulim = 1.0;
phiLim = 0.7;
km = [2 1 4];

transo = 0;
posed = zeros(1,3);

%% Simulation
phi_d = pi/4;
ad = [0; 0; -1];
zds = zf;

kp = 2; %for linear velocity
ko = 2; %for angular velocity
kh = 6; %for H gradient
kz = 7;

ky = [1.5 0; 0 1.5];
change = false;

for i = 1:samples
    
    %Manipulator control
    trans = getTransform(robot,q(:,i),end_eff_name);
    Re = trans(1:3,1:3);
    pe = trans(1:3,4);
    ze = pe(3);
    ae = Re(:,3);
    r = sqrt(ae(1:2)' * ae(1:2));
    
    ez(i) = zd(i) - ze;
    ephi(i) = phi_d - acos(ad' * ae); 
    
    fullJac = geometricJacobian(robot,q(:,i),end_eff_name);
    for j=1:ndofs
        J(:,j) = fullJac(:, dofs(j));
    end
    Jo = J(1:3,:);
    Jp = J(4:5,:);
    
    Hplot(i) = Hfunc(J);
    
    Sae = [0 -ae(3) ae(2); ae(3) 0 -ae(1); -ae(2) ae(1) 0];
    
    %Control Law #1
    Jzp = [J(6,:); (ad' * Sae * Jo)/r];
    Jzpplus = Jzp' / (Jzp * Jzp' + 0.0000 * eye(2));
    
    ud = [kz * (zd(i) - ze) + zddot(i); ko * (phi_d - atan2(r,-ae(3)))];
    qdsmall = Jzpplus * ud + (eye(4) - Jzpplus * Jzp) * ...
        (kh * gradH(dofs,J,robot,q(:,i), end_eff_name));% + 0 * kp * ep);
    

    %Mobile Robot control
    %For unicycle
    
    %To avoid atan2 dicontinuity
    td = atan2(ae(2), ae(1));
    if td - posed(3) > pi
        td = td - 2 * pi;
    elseif td - posed(3) < -pi
        td = td + 2 * pi;
    end
    
    if r < 0.01
        tdang = 0;
    else
        tdang = 1;
    end
        
    posed = [pe(1:2); tdang*td]; %desired pose (x,y,theta)
    
    ep(:,i) = -xm(1:2,i) + posed(1:2);
    emtheta(i) = xm(3,i) - posed(3);
    
    y = xm(1:2,i) + b * [cos(xm(3,i)); sin(xm(3,i))];
    yd = posed(1:2);
    yddot = Jp * qdsmall;
    
    Tthetainv = [cos(xm(3,i)) sin(xm(3,i)); -sin(xm(3,i))/b cos(xm(3,i))/b]; 
    
    ui = yddot + ky * (yd - y);
    
    res = Tthetainv * ui;
   
    %First Control Phase
    u = res(1);
    w = res(2);
    
    if norm(qdsmall) < 0.02 && norm(yd - y) < 0.001 || change == true
        change = true;
        
        %Second Control Phase
        u = 0;
        w = tdang*ad' * Sae * Sae * Jo * qdsmall / r^2 + ...
            2 * (posed(3) - xm(3,i));
    end
        
    %Update
    for j=1:ndofs
        qdot(dofs(j)) = qdsmall(j);
    end
    
    q(:,i+1) = q(:,i) + qdot*dt;
            
%     Joint limits
%     for j=1:ndofs
%        if q(dofs(j),i+1) >= maxlim(j)
%            q(dofs(j),i+1) = maxlim(j);
%        elseif q(dofs(j),i+1) <= minlim(j)
%            q(dofs(j),i+1) = minlim(j);
%        end
%     end   
    
%     u = min(max(u,-ulim),ulim);
%     w = min(max(w,-phiLim),phiLim);

    xm(1,i+1) = xm(1,i) + u * cos(xm(3,i) + w * dt / 2) * dt; %x
    xm(2,i+1) = xm(2,i) + u * sin(xm(3,i) + w * dt / 2) * dt; %y
    xm(3,i+1) = xm(3,i) + w * dt; %theta  
    
    if xm(3,i+1) > 2 * pi
        xm(3,i+1) = xm(3,i+1) - 2 * pi;
    elseif xm(3,i+1) < -2*pi
        xm(3,i+1) = xm(3,i+1) + 2 * pi;
    end  
end

%% Plotting
t = linspace(0,Time,samples);

figure;
subplot(2,2,1)
plot(t,ez);
ylabel('Z error Manipulator')
xlabel('Time')
grid on

subplot(2,2,2)
plot(t,ephi);
ylabel('phi error Manipulator')
xlabel('Time')
grid on

subplot(2,2,3:4)
plot(t,Hplot);
ylabel('H')
xlabel('Time')
grid on

figure;
subplot(2,1,1)
plot(t,ep);
ylabel('position error mobile robot')
xlabel('Time')
legend('x','y')
grid on

subplot(2,1,2)
plot(t,emtheta);
ylabel('angle error mobile robot')
xlabel('Time')
grid on

%% Gazebo Ros
step = 2; %in samples
interval = Time * step / samples;

%Ros
%rosinit
pub1 = rospublisher('/bh_j11_position_controller/command');
pub2 = rospublisher('/bh_j12_position_controller/command');
pub3 = rospublisher('/bh_j22_position_controller/command');
pub4 = rospublisher('/bh_j32_position_controller/command');
pub_arm = rospublisher('/arm_controller/command');
state_pub = rospublisher('/gazebo/set_model_state');

grip_msg = rosmessage(pub1);
pos_msg  = rosmessage(pub_arm);
state_msg = rosmessage(state_pub);

%Initial configuration
grip_msg.Data = 1.3;

pos_msg.JointNames = {'iiwa_joint_1', 'iiwa_joint_2', 'iiwa_joint_3', 'iiwa_joint_4', 'iiwa_joint_5', 'iiwa_joint_6', 'iiwa_joint_7'};
pos_msg.Points = rosmessage('trajectory_msgs/JointTrajectoryPoint');
pos_msg.Points.Positions = q(:,1);
pos_msg.Points.TimeFromStart = rosduration(1.0);

state_msg.ModelName = 'pioneer3at';
state_msg.ReferenceFrame = 'ground_plane';
state_msg.Pose.Position.X = xm(1,1);
state_msg.Pose.Position.Y = xm(2,1);
state_msg.Pose.Orientation.Z = sin(xm(3,1)/2);
state_msg.Pose.Orientation.W = cos(xm(3,1)/2);

send(pub3,grip_msg)
send(pub4,grip_msg)
send(pub_arm, pos_msg)
send(state_pub, state_msg)

pause(5)

%Sending trajectories
for i=1:samples/step
%     pos_msg.Points(i) = rosmessage('trajectory_msgs/JointTrajectoryPoint');
    pos_msg.Points.Positions = q(:,i * step);
    pos_msg.Points.TimeFromStart = rosduration(interval);
    
    state_msg.Pose.Position.X = xm(1,i * step);
    state_msg.Pose.Position.Y = xm(2,i * step);
    state_msg.Pose.Orientation.Z = sin(xm(3,i * step)/2);
    state_msg.Pose.Orientation.W = cos(xm(3,i * step)/2);
    
    send(pub_arm, pos_msg)
    send(state_pub, state_msg)
    pause(interval)
end