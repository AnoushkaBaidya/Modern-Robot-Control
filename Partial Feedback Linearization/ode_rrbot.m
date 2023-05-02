function [dX,tau1,tau2,x_d] = ode_rrbot(t,X)
    
   %initiating values of parameters 
    m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81;

    dX=zeros(4,1);  % create empty vector for [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot]
    X=num2cell(X);
    [theta1,theta2,theta1_dot,theta2_dot]=deal(X{:});

    X = zeros(4,1);
    X(1,1)= theta1;
    X(2,1)= theta2;
    X(3,1)= theta1_dot;
    X(4,1)= theta2_dot;

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

    x_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
    v_d = [q1_ddot_d;q2_ddot_d]; 

    x = [theta1;theta2;theta1_dot;theta2_dot];
    e = x - x_d;

    %initialising torque inputs 
    %Adding Gain Matrix
    k = [ 50     0    15     0;
     0     50    0     15];
    
    %tau1 =  - (k11 *theta1 +k12*theta2 + k13 *theta1_dot + k14*theta2_dot);
    %tau2 =  - (k21 *theta1 +k22*theta2 + k23 *theta1_dot + k24*theta2_dot);

    v = -k*(x-x_d) + v_d;
    
    tau = M*v + C + G;
    tau1 = tau(1);
    tau2 = tau(2);
 
    %State space representation
    %pasting values from the solve function
    dX(1)= theta1_dot;
    dX(2)= theta2_dot;
    dX(3)= (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    dX(4)= -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*tau1*cos(theta2) - 2*l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
 
 

