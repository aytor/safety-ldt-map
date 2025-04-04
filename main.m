%% AUTHOR: Aitor Gomez
%% CODE REVISED ON 05/12/2024
clear all;
clc;
%% load initial conditions and parameters
load("setup_circular.mat");
odeopts = odeset('RelTol',3e-14,'AbsTol',1e-15);
full_init_state_case = 1;

text1 = "Setting up Collision Scenario...";
fprintf(text1);

MU = 3.986*10^14; % m^3/s^2

% initial conditions
X1_0 = setup.initial_conditions(1:6);
X2_0 = setup.initial_conditions(7:12);
tspan = setup.time;

% Rotation matrices to expres covariance in Inertial frame
Rot1 = ECI_2_TNH(X1_0(1:3), X1_0(4:6));
Rot2 = ECI_2_TNH(X2_0(1:3), X2_0(4:6));

if full_init_state_case:
    % Initial full state covariance in local frame (TNH)
    Px = blkdiag(diag([1500,500,500]),diag([4.44*1e-6,0.0001,0.0001]));
    % Initial full state covariance in inertial frame (ECI)
    P1 = blkdiag(Rot1,Rot1)\Px*blkdiag(Rot1,Rot1);
    P2 = blkdiag(Rot2,Rot2)\Px*blkdiag(Rot2,Rot2);
else:
    % Initial position covariance in local frame (TNH)
    Px = blkdiag(diag([1500,1000,1000]));
    % Initial position covariance in inertial frame (ECI)
    P1 = Rot1\Px*Rot1;
    P2 = Rot2\Px*Rot2;

S0 = blkdiag(P1,P2);

dt = 1; % unscaled dt 1sec
T  = tspan(end);

% solve nominal ode
[tspan,q] = ode45(@(t,X) dynJ2(t,X), 0:dt:T, [X1_0;X2_0], odeopts);
s1 = q(:,1:6)';
s2 = q(:,7:12)';

K = length(q);
dt = (dt*1)/T;  % scaling dt down to the unit time interval [0,1] (for the free end-time OCP)

text2 = "done\n";
fprintf(text2);

%% MAP Optimization setup
import casadi.*
casadiopts = casadi.Opti();

text1 = "Setting up LDT-OCP...";
fprintf(text1);

epsilon = 1e-6;
tf = casadiopts.variable(1,1);
w1  = casadiopts.variable(3,K);
w2  = casadiopts.variable(3,K);
p1 = casadiopts.variable(6,K);
p2 = casadiopts.variable(6,K);
%dv = casadiopts.parameter(2,1);

if full_init_state_case:
    % Cost on full initial state
    obj = 0.5*epsilon*([p1(:,1);p2(:,1)]'-[s1(:,1);s2(:,1)]')*S0*([p1(:,1);p2(:,1)]-[s1(:,1);s2(:,1)]);
else:
    % Cost only on initial position
    obj = 0.5*epsilon*([p1(1:3,1);p2(1:3,1)]'-[s1(1:3,1);s2(1:3,1)]')*S0*([p1(1:3,1);p2(1:3,1)]-[s1(1:3,1);s2(1:3,1)]);
    casadiopts.subject_to(p1(4:6,1)==s1(4:6,1));
    casadiopts.subject_to(p2(4:6,1)==s2(4:6,1));

for i =1:K-1

    obj = obj + ((tf*dt)/2)*([w1(:,i);w2(:,i)]'*[w1(:,i);w2(:,i)]);

    a11 = -MU*p1(1:3,i)/norm(p1(1:3,i))^3 + w1(:,i);
    a12 = -MU*p1(1:3,i+1)/norm(p1(1:3,i+1))^3 + w1(:,i+1);

    a21 = -MU*p2(1:3,i)/norm(p2(1:3,i))^3 + w2(:,i);
    a22 = -MU*p2(1:3,i+1)/norm(p2(1:3,i+1))^3 + w2(:,i+1);

    % Integration of dynamics (leap-frog integration)
    casadiopts.subject_to( p1(1:3,i+1)==(p1(1:3,i) + (tf*dt)*p1(4:6,i) + (tf*dt)^2*0.5*a11) );
    casadiopts.subject_to( p1(4:6,i+1)==(p1(4:6,i) + (tf*dt)*0.5*(a11+a12)) );
    casadiopts.subject_to( p2(1:3,i+1)==(p2(1:3,i) + (tf*dt)*p2(4:6,i) + (tf*dt)^2*0.5*a21) );
    casadiopts.subject_to( p2(4:6,i+1)==(p2(4:6,i) + (tf*dt)*0.5*(a21+a22)) );

end

% Time constraint
casadiopts.subject_to(tf>=T-100);
casadiopts.subject_to(tf<=T+100);

% Final constraint
casadiopts.subject_to(((p2(1:3,K)-p1(1:3,K))'*(p2(1:3,K)-p1(1:3,K)))==50^2);

% Initial guesses
casadiopts.set_initial(p1,s1);
casadiopts.set_initial(p2,s2);
casadiopts.set_initial(w1,zeros(3,length(s1)));
casadiopts.set_initial(w2,zeros(3,length(s1)));
casadiopts.set_initial(tf,T-99);

casadiopts.minimize(obj);
casadiopts.solver('ipopt');
%casadiopts.ipopt.acceptable_tol =1e-7;

text2 = "done\n";
fprintf(text2);

%% solve MAP
sol = casadiopts.solve();

%% solutions
p1s = sol.value(p1);
p2s = sol.value(p2);

%% plots
figure;
hold on;

dp = 1;
p = 1e-3; % scaling factor for plots
plot3(p*s1(1,1:dp:end),p*s1(2,1:dp:end),p*s1(3,1:dp:end),'--k');
plot3(p*s2(1,1:dp:end),p*s2(2,1:dp:end),p*s2(3,1:dp:end),'--k');

plot3(p*p1s(1,1:dp:end),p*p1s(2,1:dp:end),p*p1s(3,1:dp:end),'r',LineWidth=3);
plot3(p*p2s(1,1:dp:end),p*p2s(2,1:dp:end),p*p2s(3,1:dp:end),'r',LineWidth=4);
scatter3(0,0,0,5000,[0.5,0.5,0.5],'filled')
scatter3(p*s1(1,1),p*s1(2,1),p*s1(3,1),'ok')
scatter3(p*s2(1,1),p*s2(2,1),p*s2(3,1),'ok')
scatter3(p*p1s(1,end),p*p1s(2,end),p*p1s(3,end),300,'xk')
scatter3(p*p2s(1,end),p*p2s(2,end),p*p2s(3,end),300,'xk')


xlabel('$q_x$ [m]','interpreter','latex')
ylabel('$q_y$ [m]','interpreter','latex')
zlabel('$q_z$ [m]','interpreter','latex')


axis equal;
box on;

%% zoom end point
% lim2 = 5000;
% xlim(p*[s1(1,end)-lim2, s1(1,end)+lim2])
% ylim(p*[s1(2,end)-lim2, s1(2,end)+lim2])
% zlim(p*[s1(3,end)-lim2, s1(3,end)+lim2])
%% zoom initial point s1
% lim2 = 1500;
% xlim(p*[s1(1,1)-lim2, s1(1,1)+lim2])
% ylim(p*[s1(2,1)-lim2, s1(2,1)+lim2])
% zlim(p*[s1(3,1)-lim2, s1(3,1)+lim2])
%% zoom initial point s2
% xlim(p*[s2(1,1)-lim2, s2(1,1)+lim2])
% ylim(p*[s2(2,1)-lim2, s2(2,1)+lim2])
% zlim(p*[s2(3,1)-lim2, s2(3,1)+lim2])


%% ode
function dX = dynJ2(t,X)
MU = 3.986*10^14; % m^3/s^2
Re = 6378.135;
J2 = 1082.7*1e-6;
dX = [X(4:6); ...
    -MU*X(1)/(X(1:3)'*X(1:3))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(1:3)'*X(1:3)))^2*(5*X(3)^2/(X(1:3)'*X(1:3))-1));...
    -MU*X(2)/(X(1:3)'*X(1:3))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(1:3)'*X(1:3)))^2*(5*X(3)^2/(X(1:3)'*X(1:3))-1));...
    -MU*X(3)/(X(1:3)'*X(1:3))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(1:3)'*X(1:3)))^2*(5*X(3)^2/(X(1:3)'*X(1:3))-3));
    X(10:12); ...
    -MU*X(7)/(X(7:9)'*X(7:9))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(7:9)'*X(7:9)))^2*(5*X(9)^2/(X(7:9)'*X(7:9))-1));...
    -MU*X(8)/(X(7:9)'*X(7:9))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(7:9)'*X(7:9)))^2*(5*X(9)^2/(X(7:9)'*X(7:9))-1));...
    -MU*X(9)/(X(7:9)'*X(7:9))^(3/2)*(1 - 3/2*J2*(Re/sqrt(X(7:9)'*X(7:9)))^2*(5*X(9)^2/(X(7:9)'*X(7:9))-3));];
end
