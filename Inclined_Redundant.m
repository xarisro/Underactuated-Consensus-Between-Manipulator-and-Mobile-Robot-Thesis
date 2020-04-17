%% Offline-Redundancy Algorithm and Robot Control
%Inclined_Redundant.m

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

dt = 0.01;
Time = 18;
samples = Time / dt;
t = linspace(0,Time,samples);

%% Manipulator Parameter initialization
qd = zeros(total_dofs,samples);
qd(:,1) = [pi/6; -pi/8; 0; -pi/3; 0; pi/8; pi/3]; %Initial configuration
qdot_des = zeros(total_dofs,samples);
qdot = zeros(total_dofs,1);

minlim = zeros(ndofs,1);
maxlim = zeros(ndofs,1); %Joint Limits

for i=1:ndofs
   minlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(1);
   maxlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(2);
end

J = zeros(6,ndofs);
Re = zeros(3,3);
pdot = zeros(6,1);

kp = 2; %for linear velocity
ko = 2; %for angular velocity
kh = 2; %for H gradient

%% Error initialization
ez= zeros(1,samples);
ephi = zeros(1,samples);
Hplot = zeros(1,samples);

%% Direction of apporach
trans = getTransform(robot,qd(:,1),end_eff_name);
zin = trans(3,4);
zf = 0.5;
zmax = max(zf,zin) + 0.1;
T = Time - 3;
t2 = linspace(0,T,T/dt);

w = 2 * pi / T;
zd1 = 1/2 * (1 - cos(w * t2)) * (zmax - zin) + zin;
zd2 = 1/2 * (1 + cos(w * t2)) * (zf - zmax) + zmax;

%Desired z trajectory
zd = zd1 .* (1 - heaviside(t2 - T/2)) + zd2 .* heaviside(t2 - T/2);

%Its derivative
zddot = (1/2*w*sin(w*t2)*(zmax - zin)) .* (1 - heaviside(t2 - T/2)) + ...
    (-1/2 * w * sin(w*t2)*(zf - zmax)) .* heaviside(t2 - T/2);

zd = [zd zf*ones(1,3/dt)];
zddot = [zddot zeros(1,3/dt)];

%% Simulation-Finding the Point, Redundancy Analysis
phi_d = pi/4;
ad = [0; 0; -1];
zds = zf;

for i = 1:samples
    
    %Manipulator control
    trans = getTransform(robot,qd(:,i),end_eff_name);
    Re = trans(1:3,1:3);
    p = trans(1:3,4);
    ze = p(3);
    ae = Re(:,3);
    r = sqrt(ae(1:2)' * ae(1:2));
    
    ez(i) = zds - ze;
    ephi(i) = phi_d - acos(ad' * ae);

    fullJac = geometricJacobian(robot,qd(:,i),end_eff_name);
    for j=1:ndofs
        J(:,j) = fullJac(:, dofs(j));
    end
    Jo = J(1:3,:);
    
    Hplot(i) = Hfunc(J);
    Hplot(i);
    
    % Redundant  
    Sae = [0 -ae(3) ae(2); ae(3) 0 -ae(1); -ae(2) ae(1) 0];
    
    Jnew = [J(6,:); (ad' * Sae * Jo)/r];
    Jnewplus = Jnew' / (Jnew * Jnew' + 0.0000 * eye(2));
    
    ud = [kp * (zds - ze); ko * (phi_d - atan2(r,-ae(3)))];
%     ud = [kp * (zd(i) - ze) + zddot(i); ko * (phi_d - atan2(r,-ae(3)))];
    qdsmall = Jnewplus * ud + (eye(4) - Jnewplus * Jnew) * ...
        (kh * gradH(dofs,J,robot,qd(:,i), end_eff_name));
    
    %Update
    for j=1:ndofs
        qdot(dofs(j)) = qdsmall(j);
    end
    
    qdot_des(:,i) = qdot;
    
    qd(:,i+1) = qd(:,i) + qdot*dt;
            
%     Joint limits
%     for j=1:ndofs
%        if qd(dofs(j),i+1) >= maxlim(j)
%            qd(dofs(j),i+1) = maxlim(j);
%        elseif qd(dofs(j),i+1) <= minlim(j)
%            qd(dofs(j),i+1) = minlim(j);
%        end
%     end   
end

%Final pose
qf = qd(:,samples);
qfsmall = [qf(dofs(1)); qf(dofs(2)); qf(dofs(3)); qf(dofs(4))];
tf = getTransform(robot,qf,end_eff_name);
pf = tf(1:3,4);
aef = tf(1:3,3);


%% Plotting

figure;
subplot(2,2,1)
plot(t,ephi);
ylabel('Angle error');
xlabel('Time')
grid on

subplot(2,2,2)
plot(t,ez);
ylabel('Z error')
xlabel('Time')
grid on

subplot(2,2,3:4)
plot(t,Hplot,'Linewidth',2);
ylabel('H')
xlabel('Time')
grid on

%% Unicycle Mobile Robot Parameter initialization
xm = zeros(3,samples);
xm(:,1) = [2.3; -2.3; 0.6]; %x,y,theta

ulim = 1.0;
phiLim = 0.7;
km = [0.8 2.5 3];

posed = [pf(1:2); atan2(aef(2), aef(1))]; %desired pose (x,y,theta)
transo = 0;

%% Simulation-controlling the robots

eq = zeros(1,samples);
exyz = zeros(3,samples);
ephi = zeros(1,samples);
ep = zeros(2,samples);
emtheta = zeros(1,samples);

q = zeros(total_dofs,samples);
q(:,1) = qd(:,1);
K = [2 0 0 0; 0 2 0 0; 0 0 15 0; 0 0 0 2];

for i = 1:samples
    
    %Manipulator control
    trans = getTransform(robot,q(:,i),end_eff_name);
    Re = trans(1:3,1:3);
    p = trans(1:3,4);
    ze = p(3);
    ae = Re(:,3);
    r = sqrt(ae(1:2)' * ae(1:2));
    
    exyz(:,i) = [pf(1) - p(1); pf(2) - p(2); zd(i) - ze];
    ephi(i) = phi_d - acos(ad' * ae);
    eq(i) = sqrt((qf - q(:,i))' * (qf - q(:,i))); 
    
    fullJac = geometricJacobian(robot,q(:,i),end_eff_name);
    for j=1:ndofs
        J(:,j) = fullJac(:, dofs(j));
    end
    Jo = J(1:3,:);
    
    %Inverse Jacobian - First Control Law
%     Jnew = [J(4:6,:); (ad' * Sae * Jo)/r];
%     e = [pf - p; (phi_d - atan2(r,-ae(3)))];
%     e(3) = zd(i) - ze;
%     qdsmall = Jnew' / (Jnew * Jnew' + 0.0005 * eye(4)) * ...
%         (K * e + 1*[0; 0; zddot(i); 0;]);

    %Redundant - Second Control Law
    qsmall = [q(dofs(1),i); q(dofs(2),i); q(dofs(3),i); q(dofs(4),i)];
    Jnewz = J(6,:);
    Jnewzplus = pinv(Jnewz);
    
    qdsmall = Jnewzplus * (7 * (zd(i) - ze) + zddot(i)) + ...
        (eye(4) - Jnewzplus * Jnewz) * 5 * (qfsmall - qsmall);

    %Copying A phase - Third Control Law
%     Sae = [0 -ae(3) ae(2); ae(3) 0 -ae(1); -ae(2) ae(1) 0];
%     
%     Jnew = [J(6,:); (ad' * Sae * Jo)/r];
%     Jnewplus = Jnew' / (Jnew * Jnew' + 0.0000 * eye(2));
%     
%     ud = [kp * (zd(i) - ze) + zddot(i); ko * (phi_d - atan2(r,-ae(3)))];
%     qdsmall = Jnewplus * ud + (eye(4) - Jnewplus * Jnew) * ...
%         (kh * gradH(dofs,J,robot,qd(:,i), end_eff_name));

    %Mobile Robot control
    %For unicycle
    ep(:,i) = xm(1:2,i) - posed(1:2);
    emtheta(i) = xm(3,i) - posed(3);
    
    ro = norm(ep(:,i));
    
    %To avoid atan2 dicontinuity
    transn = atan2(ep(2,i),ep(1,i));
    if transo < - pi/2 && transn > pi/2
        transn = transn - 2*pi;
    elseif transn < - pi/2 && transo > pi/2
        transn = transn + 2*pi;
    end
    transo = transn;
    
    gamma = transn - xm(3,i) + pi;
    delta = gamma + emtheta(i);
    
    u = km(1) * ro * cos(gamma);
    w = km(2) * gamma + km(1) * sin(gamma) * cos(gamma) * ...
        (gamma + km(3) * delta) / gamma;

    %Update
    for j=1:ndofs
        qdot(dofs(j)) = qdsmall(j);
    end
    
    %Trajectory following - Fourth Control Law
%     qdot = 20 * (qd(:,i) - q(:,i)) + qdot_des(:,i);
    
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
    
    if xm(3,i) > 2 * pi
        xm(3,i) = xm(3,i) - 2 * pi;
    elseif xm(3,i) < -2*pi
        xm(3,i) = xm(3,i) + 2 * pi;
    end
end

%% Plotting
t = linspace(0,Time,samples);

figure;
subplot(2,2,1)
plot(t,ephi);
ylabel('Angle error');
xlabel('Time')
grid on

subplot(2,2,2)
plot(t,exyz);
ylabel('Position error Manipulator')
xlabel('Time')
legend('x','y','z')
grid on

subplot(2,2,3:4)
plot(t,eq);
ylabel('q error Manipulator')
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