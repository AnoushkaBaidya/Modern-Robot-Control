%This is the main function
clear;clc;

%Simulation time set for 10 seconds
T=10;

%Setting initial conditions of the system 
y0 = [deg2rad(30), deg2rad(45),0,0];   % [theta_1 , theta_2, theta1_dot, thteta2_dot]

%Calling the ode45 function 
%To solve EOM
[t,y] = ode45(@ode_rrbot,[0,T],y0);

%visualize the output 
figure;
subplot(2,2,1);
plot(t,rad2deg(y(:,1)),'b','linewidth',2);
xlabel('Time in secs','FontSize',14);
ylabel('theta1 in deg','FontSize',14);
subplot(2,2,2);
plot(t,rad2deg(y(:,2)),'r','linewidth',2);
xlabel('Time in secs','FontSize',14);
ylabel('theta2 in deg','FontSize',14);
subplot(2,2,3);
plot(t,rad2deg(y(:,3)),'r','linewidth',2);
xlabel('Time in secs','FontSize',14);
ylabel('theta1dot in degrees','FontSize',14);
subplot(2,2,4);
plot(t,rad2deg(y(:,4)),'b','linewidth',2);
xlabel('Time in secs','FontSize',14);
ylabel('theta2dot in degrees','FontSize',14);

