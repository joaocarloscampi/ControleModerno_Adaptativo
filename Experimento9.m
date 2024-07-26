clc
clear 
close all

index_fig = 0;

%% Flags

saveFig = false;
%% Entrada

disp("Entrada 1: r=15");
disp("Entrada 2: r=2*sin(10*t) + 5*sin(3*t)")

entrada = input("Qual você deseja? ");
clc;

%% Experiência 9 - (d)

index_fig = 0;

% Parâmetros da planta
a = 1;
b = 5;
% Parâmetros do modelo de referência
am = 6.7;
bm = 3.3;

% Valores ideais de controle
k_star = (a+am)/b;
l_star = bm/b;

% Vetor de tempo
dt = 0.01;
t = 0:dt:7;

% Sinal de entrada
if entrada==1
    r = 15; 
    r = r*ones(1,length(t));
else
    r = 2*sin(10*t) + 5*sin(3*t);
end

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
% Sinal de entrada
if entrada==1
    gamma1 = 1;
    gamma2 = 1;
else
    gamma1 = 7;
    gamma2 = 7;
end

% Valores de x para o controlador adaptativo
x_hat_direto = zeros(1,length(t));
x0_hat = 0;
x_hat_direto(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k_direto = zeros(1,length(t));
k_direto(1) = k0;
l_direto = zeros(1,length(t));
l_direto(1) = l0;
e_direto = zeros(1,length(t));
e_direto(1) = x_hat_direto(1)-xm(1);

u_direto = zeros(1,length(t));

e = zeros(1,length(t));
e(1) = x(1)-xm(1);

u_ideal = zeros(1,length(t));

for i=2:length(t)
    % X com controlador adaptativo
    e_direto(i-1) = x_hat_direto(i-1) - xm(i-1);
    
    k_direto_dot = gamma1*e_direto(i-1)*x_hat_direto(i-1)*sign(b);
    k_direto(i) = k_direto_dot*dt + k_direto(i-1);
    
    l_direto_dot = -gamma2*e_direto(i-1)*r(i-1)*sign(b);
    l_direto(i) = l_direto_dot*dt + l_direto(i-1);
    
    u_hat_direto = -k_direto(i)*x_hat_direto(i-1) + l_direto(i)*r(i-1);
    u_direto(i) = u_hat_direto;
    
    x_hat_direto_dot = a*x_hat_direto(i-1) + b*u_hat_direto;
    x_hat_direto(i) = x_hat_direto_dot * dt + x_hat_direto(i-1);
    
    % X como controlador ideal
    e(i-1) = x(i-1) - xm(i-1);
    u = -k_star*x(i-1) + l_star*r(i-1);
    u_ideal(i) = u;
    x_dot = a*x(i-1) + b*u;
    x_i = x(i-1) + x_dot*dt;
    x(i) = x_i;
end

%% Experiencia 9 - MRAC Indireto

% Valores iniciais de adaptação
a0 = 1;
b0 = 1;
b_limite=1;
if entrada==1
    gamma1 = 1;
    gamma2 = 1;
else
    gamma1 = 7;
    gamma2 = 7;
end

% Valores de x para o controlador adaptativo
x_hat_indireto = zeros(1,length(t));
x0_hat = 0;
x_hat_indireto(1) = x0_hat;

% Vetores para progressão dos ganhos e do erro
k_indireto = zeros(1,length(t));
k_indireto(1) = k0;
l_indireto = zeros(1,length(t));
l_indireto(1) = l0;
e_indireto = zeros(1,length(t));
e_indireto(1) = x_hat_indireto(1)-xm(1);

a_indireto = zeros(1,length(t));
a_indireto(1) = a0;
b_indireto = zeros(1,length(t));
b_indireto(1) = b0;

u_indireto = zeros(1,length(t));

for i=2:length(t)
    % X com controlador adaptativo
    e_indireto(i-1) = x_hat_indireto(i-1) - xm(i-1);
    
    k_indireto_i = (am+a_indireto(i-1))/b_indireto(i-1);
    k_indireto(i) = k_indireto_i;
    
    l_indireto_i = (bm)/b_indireto(i-1);
    l_indireto(i) = l_indireto_i;
    
    u_hat_indireto = -k_indireto(i)*x_hat_indireto(i-1) + l_indireto(i)*r(i-1);
    u_indireto(i) = u_hat_indireto;
    
    a_indireto_dot = gamma1*e_indireto(i-1)*x_hat_indireto(i-1);
    a_indireto(i) = a_indireto_dot*dt + a_indireto(i-1);
    
    b_indireto_dot = gamma2*e_indireto(i-1)*u_hat_indireto;
    b_indireto(i) = b_indireto_dot*dt + b_indireto(i-1);
    
    if ~ ( (abs(b_indireto(i)) > abs(b_limite)) || ( (abs(b_indireto(i)) > abs(b_limite)) && (e_indireto(i-1)*u_hat_indireto*sgn(b) > 0) ) )
        if b_indireto(i) < 0
            b_indireto(i) = - abs(b_limite);
        else
            b_indireto(i) = abs(b_limite);
        end
    end
    
    x_hat_indireto_dot = a*x_hat_indireto(i-1) + b*u_hat_indireto;
    x_hat_indireto(i) = x_hat_indireto_dot * dt + x_hat_indireto(i-1);
end


%% Experiência 7 (a) - Plots

if entrada == 1
    titulo = "(r=15)";
else
    titulo = "(r=2sin(10t)+5sin(3t))";
end

% Plot do controlador ideal
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, x, "LineWidth", 1)
plot(t, xm, '--', "LineWidth", 2)
plot(t, e, "LineWidth", 1)
legend("Controlador Ideal", "Modelo de Referência", "Erro")
xlabel("Tempo [s]")
ylabel("Tensão")
title("Experiência 9 - Controlador ideal " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_EstadosIdeal.png')
end

% Plot do controlador direto
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, x_hat_direto, "LineWidth", 1)
plot(t, xm, '--', "LineWidth", 2)
plot(t, e_direto, "LineWidth", 1)
legend("MRAC Direto", "Modelo de Referência", "Erro")
xlabel("Tempo [s]")
ylabel("Tensão")
title("Experiência 9 - MRAC Direto " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_EstadosMRACDireto.png')
end

% Plot do controlador indireto
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, x_hat_indireto, "LineWidth", 1)
plot(t, xm, '--', "LineWidth", 2)
plot(t, e_indireto, "LineWidth", 1)
legend("MRAC Indireto", "Modelo de Referência", "Erro")
xlabel("Tempo [s]")
ylabel("Tensão")
title("Experiência 9 - MRAC Indireto " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_EstadosMRACIndireto.png')
end

% Sinal de controle
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, u_ideal, "LineWidth", 1)
plot(t, u_direto, "LineWidth", 1)
plot(t, u_indireto, "LineWidth", 1)
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]")
ylabel("Sinal u")
title("Experiência 9 - Sinal de Controle " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_SinalControle.png')
end

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, k_star_vector, '--')
plot(t, k_direto, "LineWidth", 1)
plot(t, k_indireto, "LineWidth", 1)
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]")
ylabel("Parametro")
title("Experiência 9 - Parâmetro de controle k " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_ParametroK.png')
end

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, l_star_vector, '--')
plot(t, l_direto, "LineWidth", 1)
plot(t, l_indireto, "LineWidth", 1)
legend("Ideal","MRAC Direto","MRAC Indireto")
xlabel("Tempo [s]")
ylabel("Parametro")
title("Experiência 9 - Parâmetro de controle l " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_ParametroL.png')
end

am_vetor = am*ones(1,length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, a_indireto)
plot(t, am_vetor, "LineWidth", 1)
legend("MRAC Indireto","Ideal")
xlabel("Tempo [s]")
ylabel("Parametro")
title("Experiência 9 - Parâmetro a da planta " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_ParametroA.png')
end

bm_vetor = bm*ones(1,length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, b_indireto)
plot(t, bm_vetor, "LineWidth", 1)
legend("MRAC Indireto","Ideal")
xlabel("Tempo [s]")
ylabel("Parametro")
title("Experiência 9 - Parâmetro b da planta " + titulo)
if saveFig
    saveas(gcf,'Exp9_'+string(entrada)'+'_ParametroB.png')
end