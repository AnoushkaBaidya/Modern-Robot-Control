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
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);
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
%Input Variables Array 
u1=[];
u2=[];
%Append Initial Values of the Robot at Time Zero
pos_theta1(end+1) = deg2rad(30);
pos_theta2(end+1) = deg2rad(45);
pos_theta1_dot(end+1)=0;
pos_theta2_dot(end+1)=0;
time(end+1)=0;


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

%state Matrix
X(1,1)= theta1;
X(2,1)= theta2;
X(3,1)= theta1_dot;
X(4,1)= theta2_dot;

%disp(size(X));

%Gain Matrix Values 
k11 = 23.5850;
k12 = 5.8875;
k13 = 5.1470;
k14 = 2.6108;
k21 = 5.8875;
k22 = 4.9875;
k23 = 1.5443;
k24 = 0.9770;

K =[k11,k12,k13,k14; k21,k22,k23,k24];
%Control Input to the Robot System 
u= -K*X;

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
%disp(size(pos_theta1));
%disp(size(pos_theta2));
%disp(size(pos_theta1_dot));
%disp(size(pos_theta2_dot));

%visualize the output 
figure;
subplot(2,2,1);
plot(time(1,:),rad2deg(pos_theta1(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1 in deg','FontSize',10);
subplot(2,2,2);
plot(time(1,:),rad2deg(pos_theta2(1,:)),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2 in deg','FontSize',10);
subplot(2,2,3);
plot(time(1,:),rad2deg(pos_theta1_dot(1,:)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1dot in degrees','FontSize',10);
subplot(2,2,4);
plot(time(1,:),rad2deg(pos_theta2_dot(1,:)),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2dot in degrees','FontSize',10);
figure;
subplot(2,1,1);
plot(tt(1,:),u1(1,:),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau1','FontSize',10);
subplot(2,1,2);
plot(tt(1,:),u2(1,:),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau2','FontSize',10);

