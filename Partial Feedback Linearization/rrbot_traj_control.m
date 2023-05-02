clear; close; clc;

% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);
%Time variables declaration
tic;
t = 0;
%Declaring Variables to append the position and Input Values to be plotted
%State Vectors Array to append Values 
pos_theta1=[];
pos_theta2=[];
pos_theta1_dot=[];
pos_theta2_dot=[];
%Time Variables
tt=[];
time=[];
% Torque Input Variables Array 
u1=[];
u2=[];
%Append Initial Values of the Robot at Time Zero
pos_theta1(end+1) = deg2rad(200);
pos_theta2(end+1) = deg2rad(125);
pos_theta1_dot(end+1)=0;
pos_theta2_dot(end+1)=0;
time(end+1)=0;

%Appending the desired values 
theta1_d= [];
theta2_d=[];
theta1_dot_d =[];
theta2_dot_d =[];


theta1_d(end+1) = deg2rad(200);
theta2_d(end+1) = deg2rad(125);
theta1_dot_d(end+1)=0;
theta2_dot_d(end+1)=0;


while(t < 10)
t = toc;
% read the joint states
jointData = receive(JointStates);
% inspect the "jointData" variable in MATLAB to get familiar with it structure
% implement your state feedback controller below


%Matrix X= [theta1, theta2,theta1_dot, theta2_dot]
X = zeros(4,1);
%Reading data from ROS Gazebo
theta1= jointData.Position(1,1);
theta2 = jointData.Position(2,1);
theta1_dot = jointData.Velocity(1,1);
theta2_dot= jointData.Velocity(2,1);

%storing values to plot later 
pos_theta1(end+1)= theta1;
pos_theta2(end+1)= theta2;
pos_theta1_dot(end+1)= theta1_dot;
pos_theta2_dot(end+1)= theta1_dot;
time(end+1)=t+1;
tt(end+1)=t;

%State Matrix - X
X(1,1)= theta1;
X(2,1)= theta2;
X(3,1)= theta1_dot;
X(4,1)= theta2_dot;

%disp(size(X));
% M(q),C(q,q_dot) and G(q)
M = [(9*cos(theta2))/10 + 1573/1000, (9*cos(theta2))/20 + 573/2000 ;(9*cos(theta2))/20 + 573/2000,  573/2000];
C = [-(9*theta2_dot*sin(theta2)*(2*theta1_dot + theta2_dot))/20;(9*theta1_dot^2*sin(theta2))/20];
G = [- (8829*sin(theta1 + theta2))/2000 - (28449*sin(theta1))/2000; -(8829*sin(theta1 + theta2))/2000];

%Initialising desired Trajectories 
q1_d = (pi*t^3)/500 - (3*pi*t^2)/100 + pi;
q1_dot_d = (3*pi*t^2)/500 - (3*pi*t)/50;
q1_ddot_d = (3*pi*t)/250 - (3*pi)/50;
q2_d = (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
q2_dot_d = (3*pi*t^2)/1000 - (3*pi*t)/100;
q2_ddot_d = (3*pi*t)/500 - (3*pi)/100;

theta1_d(end+1)= q1_d;
theta2_d(end+1)= q2_d;
theta1_dot_d(end+1)= q1_dot_d;
theta2_dot_d(end+1)= q2_dot_d;

 x_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
 v_d = [q1_ddot_d;q2_ddot_d]; 

 %x = [theta1;theta2;theta1_dot;theta2_dot];
 e = X - x_d;

%Gain Matrix Values 
K = [ 50 0  15  0; 0  50  0  15];


%Control Input to the Robot System 
v = -K*(X-x_d) + v_d;
u = M*v + C + G;

%Tau1 and Tau2
t1 = u(1,1);
t2 = u(2,1);

u1(end+1)=t1;
u2(end+1)=t2;

tau1.Data = t1;
tau2.Data = t2;

send(j1_effort,tau1);
send(j2_effort,tau2);
% sample the time, joint state values, and calculated torques here to be
%plotted at the end
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;
%plot the trajectories
disp(size(pos_theta1));
%disp(size(pos_theta2));
%disp(size(pos_theta1_dot));
%disp(size(pos_theta2_dot));
disp(size(theta1_d))
%visualize the output 
figure;
subplot(2,2,1);
plot(time(1,:),rad2deg(pos_theta1(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1 in deg','FontSize',10);
hold 'on';
plot(time(1,:),rad2deg(theta1_d(1,:)),'r','linewidth',2);

subplot(2,2,2);
plot(time(1,:),rad2deg(pos_theta2(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2 in deg','FontSize',10);
hold 'on';
plot(time(1,:),rad2deg(theta2_d(1,:)),'r','linewidth',2);

subplot(2,2,3);
plot(time(1,:),rad2deg(pos_theta1_dot(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1dot in degrees','FontSize',10);
hold 'on';
plot(time(1,:),rad2deg(theta1_dot_d(1,:)),'r','linewidth',2);

subplot(2,2,4);
plot(time(1,:),rad2deg(pos_theta2_dot(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2dot in degrees','FontSize',10);
hold 'on';
plot(time(1,:),rad2deg(theta2_dot_d(1,:)),'r','linewidth',2);

figure;
subplot(2,1,1);
plot(tt(1,:),u1(1,:),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau1','FontSize',10);
subplot(2,1,2);
plot(tt(1,:),u2(1,:),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau2','FontSize',10);

