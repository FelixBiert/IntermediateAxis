%% Define satellite parameters

% Inertia of the principal axes
J1 = 100;
J2 = 400;
J3 = 1600;

% angular velocity about principal axes
omega1 = 0;
omega2 = 0.01; % check for the values of Juno
omega3 = 0.0;

% Torques about principal axes
T1  = 0;
T2  = 0;
T3  = 0;

% perturbations about principal axes
epsilon1 = 0.0001;
epsilon2 = 0.000;
epsilon3 = 0.000;

tspan = [0 10000];
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