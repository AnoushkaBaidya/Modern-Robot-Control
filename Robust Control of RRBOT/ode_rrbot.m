function [dX,tau] = ode_rrbot(t,X)
    
    %Initiating Acual value parameters of th Robot
    m1 = 1; m2 = 1; l1 = 1; l2 = 1;r1 = 0.45;r2 = 0.45;I1 = 0.084; I2 = 0.084;g = 9.81;
    %Initialising Nominal Values 
    m1_hat = 0.75; m2_hat = 0.75; I1_hat = 0.063; I2_hat = 0.063;

    dX=zeros(4,1);  % create empty vector for [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot]
    X=num2cell(X);
    [theta1,theta2,theta1_dot,theta2_dot]=deal(X{:});

    X = zeros(4,1);
    X(1,1)= theta1;
    X(2,1)= theta2;
    X(3,1)= theta1_dot;
    X(4,1)= theta2_dot;

    % M(q),C(q,q_dot) and G(q) with Nominal Values-> m_hat and I_hat 
    M = [(27*cos(theta2))/40 + 4719/4000, (27*cos(theta2))/80 + 1719/8000;(27*cos(theta2))/80 + 1719/8000, 1719/8000];
    C = [-(27*theta2_dot*sin(theta2))/80, -(27*sin(theta2)*(theta1_dot + theta2_dot))/80;(27*theta1_dot*sin(theta2))/80,0];
    G = [- (26487*sin(theta1 + theta2))/8000 - (85347*sin(theta1))/8000;  -(26487*sin(theta1 + theta2))/8000];

    %Initialising desired Trajectories 
    q1_d = (pi*t^3)/500 - (3*pi*t^2)/100 + pi;
    q1_dot_d = (3*pi*t^2)/500 - (3*pi*t)/50;
    q1_ddot_d = (3*pi*t)/250 - (3*pi)/50;
    q2_d = (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
    q2_dot_d = (3*pi*t^2)/1000 - (3*pi*t)/100;
    q2_ddot_d = (3*pi*t)/500 - (3*pi)/100;

    %Desired State Trajectories
    x_d = [q1_d; q2_d; q1_dot_d;q2_dot_d];
    v_d = [q1_ddot_d;q2_ddot_d]; 

    %Current State Trajectories
    x = [theta1;theta2;theta1_dot;theta2_dot];
    e = x - x_d;

    %Initialising torque inputs for Robust Control 
    %Adding Gain Matrix
    %k = [2 0 3 0; 0 2 0 3];
    k = [12 0 7 0; 0 12 0 7];
    
    % A and B Matrix to calculate vr term
    A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    B = [0 0; 0 0; 1 0; 0 1];

    Acl = A -(B*k);
    %Calculate Vr - Robust Control Term 
    Q = eye(4).*10;
    P = lyap(Acl',Q);
    
    %Selecting a Constant Upper Bound for Uncertainity 
    %PART B with Boundary Conditions
    phi = 0.06;
    rho =15;
    %Part A without Boundary Conditions
    %rho = 7; 
    %phi = 0;
    
    %vr Conditions 
    if phi > 0
        if norm(B'*P*e) > phi
            vr = -rho*(B'*P*e)/norm(B'*P*e);
        else
            vr = -rho*(B'*P*e)/phi;
        end
    
    else
        if norm(B'*P*e) ~= 0
            vr = -rho*(B'*P*e)/norm(B'*P*e);
        else
         vr = 0;
         end
    end 

    %For part [G] remove robust control term 
    %vr = 0;
    %Virtual Input Term 
    v = -k*(e) + v_d +vr;
    %Final Control Input Term
    tau = M*v + C*[theta1_dot; theta2_dot] + G;
    tau1 = tau(1);
    tau2 = tau(2);
 
    %State space representation
    %pasting values from the solve function
    dX(1)= theta1_dot;
    dX(2)= theta2_dot;
    dX(3)= (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    dX(4)= -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*tau1*cos(theta2) - 2*l1*m2*r2*tau2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
 
 

