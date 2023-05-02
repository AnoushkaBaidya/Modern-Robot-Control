%This is the main function for PA5
clear;clc;

%RRBot system 2DOF 
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot F t g  l1 l2 r1 r2 I1 I2 tau1 tau2  'real'
syms m1 m2 'real'
%Declaring Symbolic Representation Variables
syms X1 X2 X3 X4 u lambda 'real'
syms t0 tf q0 q0_dot qf qf_dot 'real'


%Adaptive Control of RRBot 
%Generate cubic polynomial Trajectory 
%Joint1(theta1) and Joint2(theta2)

T = [1 t0 t0^2 t0^3; 0 1 2*t0 3*t0^2; 1 tf tf^2 tf^3; 0 1 2*tf 3*tf^2];
q = [q0;q0_dot;qf;qf_dot];

% a0,a1,a2,a3 
A = inv(T)*q;
A_inv = T\q; 

%Declaring the initial and final positions and velocity conditions
% Joint 1 -> [to,tf] =[0,10] [q0;q0_dot;qf;qf_dot] = [180,0,0,0]
joint_1 = subs(A,[t0,tf,q0,q0_dot,qf,qf_dot],[0,10,deg2rad(180),0,0,0]);
% Joint 1 -> [to,tf] =[0,10] [q0;q0_dot;qf;qf_dot] = [90,0,0,0]
joint_2 = subs(A,[t0,tf,q0,q0_dot,qf,qf_dot],[0,10,deg2rad(90),0,0,0]);

display(joint_1);
display(joint_2);

%joint1 q1_desired(theta1_desired) trajectory values
q1_d = joint_1(1) + joint_1(2)*t + joint_1(3)*t^2 + joint_1(4)*t^3;
q1_dot_d = jacobian(q1_d,t);
q1_ddot_d = jacobian(q1_dot_d,t);

display(q1_d);
display(q1_dot_d);
display(q1_ddot_d);

%joint2 q2_desired(theta2_desired) trajectory values
q2_d = joint_2(1) + joint_2(2)*t + joint_2(3)*t^2 + joint_2(4)*t^3;
q2_dot_d = jacobian(q2_d,t);
q2_ddot_d = jacobian(q2_dot_d,t);

display(q2_d);
display(q2_dot_d);
display(q2_ddot_d);

%q_d values and v_d values 
q_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
v_d = [q1_ddot_d;q2_ddot_d];


%A] Eqn of Motion in Linear Parametric Form 
% Regressor Matrix -> 2X5 

Y = [theta1_ddot, ...
    cos(theta2)*(2*theta1_ddot + theta2_ddot) - 2*sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
    theta2_ddot, ...
    -sin(theta1)*g, ...
    -sin(theta1 + theta2)*g; ...
    0, ...
    sin(theta2)*theta1_dot^2 + cos(theta2)*theta1_ddot, ...
    theta1_ddot + theta1_ddot, ...
    0, ...
    -sin(theta1+theta2)*g];

%Actual Values to calculate alpha actual 
%m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81;

% Alhpha Parameter Vector -> 5X1
alpha = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2;
        m2*l1*r2
        m2*r2^2 + I2
        m1*r1 + m2*l1
        m2*r2];

%Eqn of Motion in Linear parametric Form 
tau_LinearPara= Y*alpha;
display(tau_LinearPara)

alpha = subs(alpha,[m1,m2,l1,l2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);
display(alpha);

%Comparing with Manipulator Equation form 
a = I1 + I2 + m1*r1^2 + m2*(l1^2 + r2^2);
b = m2*l1*r2;
d = I2 + m2*r2^2;

Mmat= [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
Cmat= [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
Gmat= [-m1*g*r1*sin(theta1)-m2*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2*g*r2*sin(theta1+theta2)];

%display(Mmat);
%display(Cmat);
%display(Gmat);



%Control Input in Manipulator Form 
tau_Manipulator = Mmat*[theta1_ddot;theta2_ddot] + Cmat*[theta1_dot;theta2_dot] + Gmat;
display(tau_Manipulator)

%Alpha hat values -> Nominal Unknown Parameter Values 
%Initiating Nominal value parameters of th Robot
m1 = 0.75; m2 = 0.75; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.063; I2 = 0.063;g = 9.81;
alpha_hat = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2;
        m2*l1*r2
        m2*r2^2 + I2
        m1*r1 + m2*l1
        m2*r2];

display(alpha_hat)

%using State Feeedback Control Design to determine control Gains Kp and Kd
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

%Step1 :Check the controllability of the system 
Rank_Matrix = rank(ctrb(A,B));
%Rank = 4 which is greater than n thus Controllable System

% State Feedback Design for Eigen Values {-1,-1,-2,-2}
lambda = [-2, -2, -3, -3];
%lambda =[-3,-3,-4,-4];
%kn = [Kp,Kd]
Kn = place(A,B,lambda);
Kp = [Kn(:,1),Kn(:,2)];
Kd = [Kn(:,3),Kn(:,4)];

sz_K = size(Kn); % kn = 2X4

%Acl -> A closed Loop Solution -> 2nX2n
Acl = [zeros(2),eye(2);-Kp,-Kd]; 
display(Acl);

%Tuning Variables for Tuning the Adaptive COntroller
Q = eye(4);
P = lyap(Acl', Q);
%Checking If P is Positive Definite(PD) matrix
Eig_P = eig(P);
%Display Eigenvalues of P Matrix
display(Eig_P);

gamma = eye(5).*0.1;
T=10;

%Calling the ode45 function
%Simulate Using ode45
x0 = [deg2rad(200); deg2rad(125);0;0;alpha_hat(1);alpha_hat(2);alpha_hat(3);alpha_hat(4);alpha_hat(5)];
[t, y] = ode45(@(t,X) ode_rrbot(t,X,Kn,P,gamma,B), [0,T], x0, [0;0]);


u = zeros(length(t),2);
theta1_d_list = [];
theta2_d_list =[];
theta1_dot_d_list =[];
theta2_dot_d_list =[];

%Reconstruct Control Input 
for i = 1:length(t)
    [~,u(i,:)] = ode_rrbot(t(i,:),y(i,:),Kn,P,gamma,B);
    theta1_d = (pi*t(i)^3)/500 - (3*pi*t(i)^2)/100 + pi;
    theta2_d = (pi*t(i)^3)/1000 - (3*pi*t(i)^2)/200 + pi/2;
    theta1_dot_d = (3*pi*t(i)^2)/500 - (3*pi*t(i))/50;
    theta2_dot_d = (3*pi*t(i)^2)/1000 - (3*pi*t(i))/100;
    theta1_ddot_d = (3*pi*t(i))/250 - (3*pi)/50;
    theta2_ddot_d = (3*pi*t(i))/500 - (3*pi)/100;

    theta1_d_list = [theta1_d_list,theta1_d];
    theta2_d_list = [theta2_d_list,theta2_d];
    theta1_dot_d_list = [theta1_dot_d_list,theta1_dot_d];
    theta2_dot_d_list = [theta2_dot_d_list,theta2_dot_d];


end 

%Plot the Data 
%visualize the output 
figure;
text(0.5, 100, 'Blue: State Trajectory and Red: Desired State Trajectory', 'HorizontalAlignment', 'center', 'FontSize', 12);
subplot(2,2,1);
plot(t,rad2deg(y(:,1)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1 in deg','FontSize',10);
hold 'on';
plot(t,rad2deg(theta1_d_list(1,:)),'r','linewidth',2);

subplot(2,2,2);
plot(t,rad2deg(y(:,2)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2 in deg','FontSize',10);
hold 'on';
plot(t,rad2deg(theta2_d_list(1,:)),'r','linewidth',2);

subplot(2,2,3);
plot(t,rad2deg(y(:,3)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta1dot in degrees','FontSize',10);
hold 'on';
plot(t,rad2deg(theta1_dot_d_list(1,:)),'r','linewidth',2);

subplot(2,2,4);
plot(t,rad2deg(y(:,4)),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('theta2dot in degrees','FontSize',10);
hold 'on';
plot(t,rad2deg(theta2_dot_d_list(1,:)),'r','linewidth',2);

figure;
subplot(2,1,1);
plot(t,u(:,1),'b','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau1','FontSize',10);
subplot(2,1,2);
plot(t,u(:,2),'r','linewidth',2);
xlabel('Time in secs','FontSize',10);
ylabel('tau2','FontSize',10);

figure;
colors = lines(5); % choose 5 distinct colors for each alpha hat and alpha pair
plot(t,y(:,5:9),'linewidth',2);
xlabel('time (sec)');
ylabel('alphahat');
hold on;
for i = 1:5
    plot(t,alpha(i)*ones(1,length(t)),'--','linewidth',2,'color',colors(i,:));
    % set the color of each alpha line to the same color as its corresponding alpha hat line
end
legend('alphahat 1', 'alphahat 2', 'alphahat 3', 'alphahat 4', 'alphahat 5', 'alpha 1' , 'alpha 2' ,'alpha 3' ,'alpha 4','alpha 5');
