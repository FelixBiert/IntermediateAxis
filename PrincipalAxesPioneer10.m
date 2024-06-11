%% ---------------------------- Description ------------------------------- 
% Simple example to show the effects of spinning a spacecraft about its
% principal axes
% Spacecraft: Pioneer 10
% Spin rate: 4.8 rpm
% Spin axis: major axis during mission
% Detailed spacecraft data provided in:
% Pioneer F/G: Spacecraft Operational Characteristics (NASA/TP-2019-220324)

%% Define satellite parameters

% Inertia of the principal axes
J1 = 284.5*1.36; % conversion from slug ft^2 to kg m^2
J2 = 180.5*1.36; % conversion from slug ft^2 to kg m^2
J3 = 433.9*1.36; % conversion from slug ft^2 to kg m^2

% angular velocity about principal axes
omega1 = 0.5027;
omega2 = 0.0; % [rad/s]
omega3 = 0.0;

% Torques about principal axes
T1  = 0;
T2  = 0;
T3  = 0;

% perturbations about principal axes
epsilon1 = 0.001;
epsilon2 = 0.001;
epsilon3 = 0.001;

tspan = [0 600];
x0 = [omega1+epsilon1, omega2+epsilon2, omega3+epsilon3];
[t,x] = ode45(@(t,x) odefcn(t,x,J1,J2, J3, T1, T2, T3), tspan, x0);

%% plot the results
figure('Name', 'Group plot')
subplot(3,1,1)
plot(t, x(:, 1))
xlabel('Time [s]')
ylabel('\omega_1 [rad/s]')
ylim([-0.03, 0.03])

subplot(3,1,2)
plot(t, x(:, 2))
xlabel('Time [s]')
ylabel('\omega_2 [rad/s]')
ylim([-0.03, 0.03])

subplot(3,1,3)
plot(t, x(:, 3))
xlabel('Time [s]')
ylabel('\omega_3 [rad/s]')
ylim([-0.03, 0.03])

figure('Name', 'All in one')

plot(t, x(:, 1))
hold on
plot(t, x(:, 2))
plot(t, x(:, 3))

xlabel('Time [s]')

%% Nonlinaer Equations of Motion
function dxdt = odefcn(t,x,J1,J2, J3, T1, T2, T3)
  
dxdt = zeros(3,1);

dxdt(1) = T1/J1-((J3-J2)/J1)*x(2)*x(3);
dxdt(2) = T2/J2-((J1-J3)/J2)*x(1)*x(3);
dxdt(3) = T3/J3-((J2-J1)/J3)*x(1)*x(2);

end