%% Offline-Langrange Algorithm
%Langrange.m

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

Nmax = 3000; %Algorithm Iterations

%% Manipulator Parameter initialization
q = [pi/3; -pi/8; 0; -pi/2; 0; pi/8; pi/3]; %Initial configuration
qdot = zeros(total_dofs,1);

minlim = zeros(ndofs,1);
maxlim = zeros(ndofs,1); %Joint Limits

for i=1:ndofs
   minlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(1);
   maxlim(i) = robot.Bodies{dofs(i)+1}.Joint.PositionLimits(2);
end

J = zeros(6,ndofs);
Re = zeros(3,3);

%% Errors initialization
ez= zeros(1,Nmax);
etheta = zeros(1,Nmax);
Hplot = zeros(1,Nmax);

%% Simulation-Finding the Point, Lagrange multipliers
a = 0.03;
theta_d = pi / 4;
zd = 0.4;
ad = [0; 0; -1];
lambda = [0.1; -0.1];

kz = 1;
ka = 1;

if theta_d < 0.1
    ka = 7;
end

for i = 1:Nmax
    
    %Manipulator kinematics
    trans = getTransform(robot,q,end_eff_name);
    Re = trans(1:3,1:3);
    p = trans(1:3,4);
    ze = p(3);
    ae = Re(:,3);
    
    ez(i) = zd - ze;
    etheta(i) = theta_d - acos(ad' * ae);

    fullJac = geometricJacobian(robot,q,end_eff_name);
    for j=1:ndofs
        J(:,j) = fullJac(:, dofs(j));
    end
    Jo = J(1:3,:);
    
    Hplot(i) = Hfunc(J);
%     Hplot(i);
    
    % Lagrangian  
    g1 = kz * (ze - zd);
    g2 = ka * (ad' * ae - cos(theta_d));
    
    Sae = [0 -ae(3) ae(2); ae(3) 0 -ae(1); -ae(2) ae(1) 0];
    gradg1 = kz * J(6,:)';
    gradg2 = ka * (-ad' * Sae * Jo)';

    Lx = -gradH(dofs,J,robot,q, end_eff_name) + [gradg1 gradg2] * lambda;
    Ll = [g1; g2];
   
    
    %Update
    for j=1:ndofs
        q(dofs(j)) = q(dofs(j)) - a * Lx(j);
    end
    
    lambda = lambda + a * Ll;
        
%     Joint limits
%     for j=1:ndofs
%        if q(dofs(j)) >= maxlim(j)
%            q(dofs(j)) = maxlim(j);
%        elseif q(dofs(j)) <= minlim(j)
%            q(dofs(j)) = minlim(j);
%        end
%     end   
end

%% Plotting

figure;
subplot(2,2,1)
plot(etheta);
subplot
ylabel('Angle error');
xlabel('Number of Iterations')
grid on

subplot(2,2,2)
plot(ez);
ylabel('Z error')
xlabel('Number of Iterations')
grid on

subplot(2,2,3:4)
plot(Hplot);
ylabel('H')
xlabel('Number of Iterations')
grid on
