clc
clear 
close all

index_fig = 0;

%% Experiência 7 - (a)

clear;
clc;
close all;
index_fig = 0;

% Parâmetros da planta
a = 2;
b = 5;
% Parâmetros do modelo de referência
am = 10;
bm = 1;

% Valores ideais de controle
k_star = (a+am)/b;
l_star = bm/b;

% Vetor de tempo
dt = 0.01;
t = 0:dt:10;

% Sinal de entrada
r = 5*cos(t) + 2*sin(pi*t);

% Valores de x para o controlador ideal
x0 = 0;
x = zeros(1,length(t));
x(1) = x0;

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);
xm = lsim(G, r, t);

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
gamma1 = 50;
gamma2 = 0.5;

% Valores de x para o controlador adaptativo
x_hat = zeros(1,length(t));
x0_hat = 0;
x_hat(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(1,length(t));
k(1) = k0;
l = zeros(1,length(t));
l(1) = l0;
e = zeros(1,length(t));
e(1) = x_hat(1)-xm(1);

for i=2:length(t)
    % X com controlador adaptativo
    e(i-1) = x_hat(i-1) - xm(i-1);
    
    k_dot = gamma1*e(i-1)*x_hat(i-1)*sign(b);
    k(i) = k_dot*dt + k(i-1);
    
    l_dot = -gamma2*e(i-1)*r(i-1)*sign(b);
    l(i) = l_dot*dt + l(i-1);
    
    u_hat = -k(i)*x_hat(i-1) + l(i)*r(i-1);
    x_hat_dot = a*x_hat(i-1) + b*u_hat;
    x_hat(i) = x_hat_dot * dt + x_hat(i-1);
    
    % X como controlador ideal
    u = -k_star*x(i-1) + l_star*r(i-1);
    x_dot = a*x(i-1) + b*u;
    x_i = x(i-1) + x_dot*dt;
    x(i) = x_i;
end

%% Experiência 7 (a) - Plots

figure(1)
hold on
grid on
plot(t, x, "LineWidth", 1)
plot(t, xm, '--', "LineWidth", 2)
plot(t, x_hat, "LineWidth", 1)
legend("Controlador Ideal", "Modelo de Referência", "Adaptativo")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(a) - Estado x")


figure(2)
hold on
grid on
plot(t,e)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Erro de adaptação")

k_star_vector = k_star*ones(1, length(t));
figure(3)
hold on
grid on
plot(t,k)
plot(t, k_star_vector, '--')
legend("Adaptativo", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Parâmetro de controle k")
ylim([0, 9])

l_star_vector = l_star*ones(1, length(t));
figure(4)
hold on
grid on
plot(t,l)
plot(t, l_star_vector, '--')
legend("Adaptativo", "Ideal")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(a) - Parâmetro de controle l")
ylim([0, 1.5])