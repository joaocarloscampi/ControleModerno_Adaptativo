close all;
clc;
clear;

%% Flags

saveFig = false;

%% Experiência 7 (d)

index_fig = 0;

% Parametros da planta
am = 10;
bm = 10;

a = -5;
b =  3;

% Valores ideais de controle
k_star = (am-a)/b;
l_star = bm/b;

% Intervalo de tempo
dt = 0.01;
t = 0:dt:10;

% Valores de referência
r = ones(1,length(t));

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);

ym = lsim(G, r, t);

% Vetores para adaptação
u = zeros(1,length(t));
uf = zeros(1,length(t));
theta = zeros(2,length(t));
phi = zeros(2,length(t));
xi = zeros(1,length(t));
b_hat = zeros(1,length(t));
ms_2 = zeros(1,length(t));
epslon = zeros(1,length(t));
yp = zeros(1,length(t));

% Valores iniciais
Gamma = [200, 50;
         50, 50];
     
gamma = 10;

yp0 = 0;
k0 = 1;
l0 = 1;
uf0 = 0;
xi_0 = [k0, l0]*[yp0; r(1)] + uf0;
b_hat_0 = 0;

yp(1) = yp0;
theta(:,1) = [k0;l0];
uf(1) = uf0;
xi(1) = xi_0;
b_hat(1) = b_hat_0;
phi(:,1) = [yp0; -r(1)];


for i=2:length(t)
    uf_dot = -am*uf(i-1) + u(i-1);
    uf(i) = uf_dot*dt + uf(i-1);
    
    phi_dot_1 = -am*phi(1,i-1) + yp(i-1);
    phi_dot_2 = -am*phi(2,i-1) - r(i-1);
    phi(1,i) = phi_dot_1*dt + phi(1,i-1);
    phi(2,i) = phi_dot_2*dt + phi(2,i-1);
    
    xi(i) = theta(:,i-1)'*phi(:,i) + uf(i); % corrigir theta aqui
    
    b_hat_dot = gamma*epslon(i-1)*xi(i-1);
    b_hat(i) = b_hat_dot*dt + b_hat(i-1);
    
    e(i) = yp(i-1) - ym(i);
    ms_2(i) = 1 + phi(:,i)'*phi(:,i) + uf(i)^2;
    epslon(i) = (e(i)-b_hat(i)*xi(i))/ms_2(i);
    
    theta_dot = Gamma*epslon(i)*phi(:,i);
    theta(:,i) = theta_dot*dt + theta(:,i-1);
    
    k_i = theta(1,i);
    l_i = theta(2,i);
    u(i) = -k_i*yp(i-1) + l_i*r(i-1);
    
    yp_dot = -a*yp(i-1) + b*u(i);
    yp(i) = yp_dot*dt + yp(i-1);
end

%% Plots

% Plot de x

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, yp)
plot(t, ym)
legend("Estimado", "Modelo de Referência")
xlabel("Tempo [s]")
ylabel("Amplitude")
title("Experiência 7(d) - Saída y")
if saveFig
    saveas(gcf,'Ex7D_SaidaY.png')
end

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, e)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(d) - Erro")
if saveFig
    saveas(gcf,'Ex7D_Erro.png')
end

% Plot de k

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, theta(1,:))
plot(t, k_star_vector)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(d) - Parametro k")
legend("Estimado", "Ideal")
if saveFig
    saveas(gcf,'Ex7D_ParametroK.png')
end

% Plot de k

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, theta(2,:))
plot(t, l_star_vector)
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(d) - Parametro l")
legend("Estimado", "Ideal")
if saveFig
    saveas(gcf,'Ex7D_ParametroL.png')
end
