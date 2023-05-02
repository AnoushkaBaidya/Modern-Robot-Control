%This is the main function
clear;clc;

%RRBot system 2DOF 
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot F t g  l1 l2 r1 r2 I1 I2 tau1 tau2  'real'
syms m1 m2 'real'
%Declaring Symbolic Representation Variables
syms X1 X2 X3 X4 u lambda 'real'

%Calculate the Langrangian Equation L=K-P
K = (1/2)*m1*(r1^2)*(theta1_dot^2) + (1/2)*I1*(theta1_dot^2) + (1/2)*m2*(l1^2)*(theta1_dot^2) + (1/2)*m2*(r2^2)*((theta2_dot+theta1_dot)^2) + m2*l1*r2*theta1_dot*(theta1_dot+theta2_dot)*cos(theta2) +(1/2)*I2*((theta2_dot+theta1_dot)^2);
%display(K)

P = m1*g*r1*cos(theta1) + m2*g*l1*cos(theta1) + m2*g*r2*cos(theta1 + theta2);
%display(P)

%Lagrangian Equation 
L=K-P;
%display(L)

%Solving for Step4 
dL_dtheta1 = jacobian(L,theta1);
dL_dtheta2 = jacobian(L,theta2);

dL_dtheta1_dot = jacobian(L, theta1_dot);
dL_dtheta2_dot = jacobian(L, theta2_dot);

ddt_dL_dtheta1_dot = jacobian(dL_dtheta1_dot,[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];
ddt_dL_dtheta2_dot = jacobian(dL_dtheta2_dot,[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];
%display(ddt_dL_dtheta1_dot);
%display(ddt_dL_dtheta2_dot);

%Equations of Motion 
eqn1=ddt_dL_dtheta1_dot - dL_dtheta1 -tau1;
eqn2=ddt_dL_dtheta2_dot - dL_dtheta2 -tau2;
%display(eq1)
%display(eq2)

sol = solve([eqn1==0, eqn2==0], [theta1_ddot, theta2_ddot]);

%display(sol.theta1_ddot)
%display(sol.theta2_ddot)

%A] Finding the Equilibrium Points for the system
eqn = [eqn1;eqn2];
EOM = subs(eqn, [theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,tau1,tau2],[0,0,0,0,0,0]);

sol = solve(EOM == 0,[theta1,theta2]);

display(sol.theta1)
display(sol.theta2)

%Defining Physical Properties of the system 
 m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81; 

%State Space Representation
%State Vector and Generalised Input Vector 
X = [X1,X2,X3,X4];
u = [tau1;tau2];

X1_dot = X3;
X2_dot = X4;
X3_dot = (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
X4_dot = -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*tau1*cos(theta2) - 2*l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

X3_dot= subs(X3_dot,[theta1,theta2,theta1_dot,theta2_dot],[X1,X2,0,0]);
X4_dot= subs(X4_dot,[theta1,theta2,theta1_dot,theta2_dot],[X1,X2,0,0]);

X_dot = [X1_dot;X2_dot;X3_dot;X4_dot];

% B] Jacobian Linearization Around Equilibrium Point 
% A = State Matrix
A = jacobian(X_dot,X);
%B = Input Matrix
B = jacobian(X_dot,u);

%Checking the size of A and B matrix 
sz_B = size(B);
sz_A = size(A);

%Linearize around upward equilibrium point [0,0,0,0]
A1 = subs(A, [X1,X2,X3,X4],[0,0,0,0]);
B1 = subs(B, [X1,X2,X3,X4],[0,0,0,0]);

A1 = double(A1);
B1 = double(B1);

%C] Checking for Stability of the system 
eig_A1 = eig(A1);
disp(eig_A1);

%Linearize around upward equilibrium point [pi,0,0,0]
A2 = subs(A, [X1,X2,X3,X4],[pi,0,0,0]);
B2 = subs(B, [X1,X2,X3,X4],[pi,0,0,0]);

A2 = double(A2);
B2 = double(B2);

% Checking for Stability of the system 
eig_A2 = eig(A2);
disp(eig_A2);

%Linearize around upward equilibrium point [0,pi,0,0]
A3 = subs(A, [X1,X2,X3,X4],[0,pi,0,0]);
B3 = subs(B, [X1,X2,X3,X4],[0,pi,0,0]);

A3 = double(A3);
B3 = double(B3);

% Checking for Stability of the system 
eig_A3 = eig(A3);
disp(eig_A3);

%Design a controller to control the system 
%D] Step1: Check Controllability of the system
%Check Rank of C
Rank_Matrix = rank(ctrb(A1,B1))
%Rank = 4 which is greater than n thus Controllable System 

%E] State Feedback Design for Eigen Values of my choice 
lambda = [-1,-2,-1+1i, -1-1i];

Kn = place(A1,B1,lambda);

sz_K = size(Kn);

%Step3 of Designing Controller 


%Simulation time set for 10 seconds
T=10;

%Setting initial conditions of the system 
y0 = [deg2rad(30), deg2rad(45),0,0];   % [theta_1 , theta_2, theta1_dot, thteta2_dot]

%Calling the ode45 function 
%Simulate Using ode45
[t,y] = ode45(@ode_rrbot,[0,T],y0);

%Reconstruct Control Input 
%Adding Gain Matrix
k11 = 23.5850;
k12 = 5.8875;
k13 = 5.1470;
k14 = 2.6108;
k21 = 5.8875;
k22 = 4.9875;
k23 = 1.5443;
k24 = 0.9770;

K =[k11,k12,k13,k14; k21,k22,k23,k24];
%disp(K)

%Control Input 
u = -K*y';

t1=u(1,:);
t2=u(2,:);


%visualize the output 
figure;
subplot(2,3,1);
plot(t,rad2deg(y(:,1)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1 in deg','FontSize',10);
subplot(2,3,2);
plot(t,rad2deg(y(:,2)),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2 in deg','FontSize',10);
subplot(2,3,3);
plot(t,rad2deg(y(:,3)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1dot in degrees','FontSize',10);
subplot(2,3,4);
plot(t,rad2deg(y(:,4)),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2dot in degrees','FontSize',10);
subplot(2,3,5);
plot(t,u(1,:),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau1','FontSize',10);
subplot(2,3,6);
plot(t,u(2,:),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau2','FontSize',10);






