% clear all
clf

a = 1; % The last digit of the day you were born (13.4.1990 -> a = 3)
        %set 1 if the date ends with zero
b = 1; % The last digit of the month you were born (13.4.1990 -> b = 4)
        %set 1 if the month ends with zero

c = 4*((a+b-2)/16+1);
num = [1 3/2*c]; % numerator of the transfer function of the system
den = conv([1 c], [1 -1]); % denominator of the transfer function of the system
G = tf(num, den)
[A, B, C, D] = tf2ss(num, den);
x0 = [0 0.1]'; %initial conditions

Tend = 5; % Final time of the simulation
Ts = 0.05; % Control period

u_sat = 1.5; %saturation limits

%% Simulation
t = 0:Ts:Tend;
u = zeros(numel(t), 1);
y = zeros(numel(t), 1);
x_sys = zeros(numel(t)+1, size(A,1));
x_sys(1,:) = x0;

r = t>=1; %reference signal - step at time 1


e = zeros(numel(t), 1);
ssum = zeros(numel(t), 1);

kp = 9;
ks = 2;

for k = 1:numel(t)
    y(k) = C*x_sys(k,:)'; % output of the system

    % == Your code goes here ==
    e(k) = r(k) - y(k); % PS controller input
    ssum(k + 1) = ssum(k) + e(k);
    
    u(k) = kp*e(k) + ks*ssum(k+1); 
   
    if abs(u(k)) >= u_sat % anti-windup
       u(k) = kp*e(k);
       ssum(k+1) = 0;
    end
   
    % Saturation of the input
    u(k) = max(min(u(k), u_sat), -u_sat);
    
    % Simulation from time t(k) to time t(k+1)
    [~, x_tmp] = ode45(@(t, x) A*x + B*u(k), [0; Ts], x_sys(k,:)');
    x_sys(k+1,:) = x_tmp(end,:); %new state at time k+1
end

%% Visualization
figure(1)
clf

subplot(211)
stairs(t, r, 'k--')
hold on
plot(t, y)
hold off
grid on
xlabel('Time [s]')
ylabel('Output y(t)')
title('Response of the system')

subplot(212)
stairs(t, u)
hold on
plot([0 Tend], u_sat*[1 1], 'k-.')
plot([0 Tend], -u_sat*[1 1], 'k-.')
hold off
grid on
ylim(1.2*[-u_sat u_sat])
xlabel('Time [s]')
ylabel('Input u(t)')

% checks if last 10 samples have (roughly) the same value
if mean(abs(diff(y(end-10 : end)))) > 0.05
    error('System is not stabilized.')
end

title('Output of the regulator')

if sum(abs(y-1) <= 0.05)
    I = find(abs(y-1) > 0.05);
    settling_time = t(I(end)+1) - 1;
else
    settling_time = Inf;
end

fprintf('Overshoot: %.1f %%\n', (max(y)-1)*100)
fprintf('Settling time: %.1f s\n', settling_time)



