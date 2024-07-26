clc
clear 
close all

index_fig = 0;

%% Experiência 7 - (b)

close all;
clc;
clear;
index_fig = 0;

a = -1;
b = 12;
am = 2;
bm = 2;

% Valores ideais de controle
k_star = (a+am)/b;
l_star = bm/b;

% Intervalo de tempo
dt = 0.01;
t = 0:dt:10;

% Valores de referência
r1 = 5*ones(1,length(t));
r2 = sin(2*t);
r3 = 1./(1+t);

r = [r1; r2; r3];

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);

xm_1 = lsim(G, r1, t);
xm_2 = lsim(G, r2, t);
xm_3 = lsim(G, r3, t);

xm = [xm_1 xm_2 xm_3]';

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
gamma1 = 1;
gamma2 = 20;

% Valores de x para o controlador adaptativo
x_hat = zeros(3,length(t));
x0_hat = 0;
x_hat(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(3,length(t));
k(1) = k0;
l = zeros(3,length(t));
l(1) = l0;
e = zeros(3,length(t));
e(1) = x_hat(1)-xm(1);

for j=[1,2,3]
    for i=2:length(t)
        % X com controlador adaptativo
        e(j,i-1) = x_hat(j,i-1) - xm(j,i-1);

        k_dot = gamma1*e(j,i-1)*x_hat(j,i-1)*sign(b);
        k(j,i) = k_dot*dt + k(j,i-1);

        l_dot = -gamma2*e(j,i-1)*r(j,i-1)*sign(b);
        l(j,i) = l_dot*dt + l(j,i-1);

        u_hat = -k(j,i)*x_hat(j,i-1) + l(j,i)*r(j,i-1);
        x_hat_dot = a*x_hat(j,i-1) + b*u_hat;
        x_hat(j,i) = x_hat_dot * dt + x_hat(j,i-1);
    end
end

%% Experiência 7 (b) - Plots

% Plot de x

index_fig = index_fig + 1;
figure(index_fig)
subplot(3,1,1)
hold on
grid on
plot(t, xm(1,:), '--', "LineWidth", 2)
plot(t, x_hat(1,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b): r=5")
subplot(3,1,2)
hold on
grid on
plot(t, xm(2,:), '--', "LineWidth", 2)
plot(t, x_hat(2,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b) - sin(2t)")
subplot(3,1,3)
hold on
grid on
plot(t, xm(3,:), '--', "LineWidth", 2)
plot(t, x_hat(3,:), "LineWidth", 1)
legend("Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(b)  r=1/(1+t)")

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
subplot(3,1,1)
hold on
grid on
plot(t,e(1,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b): r=5")
subplot(3,1,2)
hold on
grid on
plot(t,e(2,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - sin(2t)")
subplot(3,1,3)
hold on
grid on
plot(t,e(3,:))
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b)  r=1/(1+t)")

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,k)
plot(t, k_star_vector, '--')
legend("r=5", "r=sin(2t)", "r=1/(1+t)", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - Parâmetro de controle k")

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l)
plot(t, l_star_vector, '--')
legend("r=5", "r=sin(2t)", "r=1/(1+t)", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(b) - Parâmetro de controle L")