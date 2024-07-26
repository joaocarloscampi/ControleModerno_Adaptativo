close all;
clc;
clear;

%% Flags

saveFig = false;

%% Experiência 7 (c) - Direto (i) 

index_fig = 0;

a = 0.5;
b = 1.5;
d = 10;

am = 0.5;
bm = 0.5;

% Valores ideais de controle
k_star = (am-a)/b;
l_star = bm/b;
delta_star = -d*a/b;

% Intervalo de tempo
dt = 0.01;
t = 0:dt:10;

% Valores de referência
Vs = 55*ones(1,length(t));

% Função de transferência modelo de referência
s = tf('s');
G = bm/(s+am);

Vm = lsim(G, Vs, t);

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
delta0 = 1;

gamma1 = 2;
gamma2 = 1;
gamma3 = 1;

% Valores de V para o controlador adaptativo
V_hat = zeros(1,length(t));
V0_hat = 0;
V_hat(1) = V0_hat;

% Vetores para progressão dos ganhos e do erro
k = zeros(1,length(t));
k(1) = k0;
l = zeros(1,length(t));
l(1) = l0;
delta = zeros(1,length(t));
delta(1) = delta0;

e = zeros(1,length(t));
e(1) = V_hat(1)-Vm(1);

for i=2:length(t)
    % X com controlador adaptativo
    e(i-1) = V_hat(i-1) - Vm(i-1);

    k_dot = gamma1*e(i-1)*V_hat(i-1)*sign(b);
    k(i) = k_dot*dt + k(i-1);

    l_dot = -gamma2*e(i-1)*Vs(i-1)*sign(b);
    l(i) = l_dot*dt + l(i-1);
    
    delta_dot = -gamma3*e(i-1)*sign(b);
    delta(i) = delta_dot*dt + delta(i-1);

    u_hat = -k(i)*V_hat(i-1) + l(i)*Vs(i-1) + delta(i);
    V_hat_dot = - a*V_hat(i-1) + b*u_hat + d*a;
    V_hat(i) = V_hat_dot * dt + V_hat(i-1);
end


%% Experiência 7 (c) - Direto (ii)
%  RODAR SEÇÃO ANTERIOR ANTES!

% Valores iniciais de adaptação
k0 = 1;
l0 = 1;
delta0 = 1;

gamma1 = 2;
gamma2 = 1;
gamma3 = 1;

V0_hat = 0;

% Valores de V para o controlador adaptativo
V_hat_ii = zeros(1,length(t));
V_hat_ii(1) = V0_hat;

% Vetores para progressão dos ganhos e do erro
k_ii = zeros(1,length(t));
k_ii(1) = k0;
l_ii = zeros(1,length(t));
l_ii(1) = l0;
delta_ii = zeros(1,length(t));
delta_ii(1) = delta0;

e_ii = zeros(1,length(t));
e_ii(1) = V_hat_ii(1)-Vm(1);

for i=2:length(t)
    % X com controlador adaptativo
    e_ii(i-1) = V_hat_ii(i-1) - Vm(i-1);

    k_ii_dot = gamma1*e_ii(i-1)*V_hat_ii(i-1)*sign(b);
    k_ii(i) = k_ii_dot*dt + k_ii(i-1);

    l_ii_dot = -gamma2*e_ii(i-1)*Vs(i-1)*sign(b);
    l_ii(i) = l_ii_dot*dt + l_ii(i-1);
    
    delta_ii_dot = -gamma3*e_ii(i-1)*sign(b);
    delta_ii(i) = delta_ii_dot*dt + delta_ii(i-1);

    u_hat = -k_ii(i)*V_hat_ii(i-1) + l_ii(i)*Vs(i-1) + delta_ii(i);
    a1 = a + 0.04/(1+V_hat_ii(i-1));
    d1 = 0.2 + sin(0.02*t(i-1));
    V_hat_dot = - a1*V_hat_ii(i-1) + b*u_hat + d1*a1;
    V_hat_ii(i) = V_hat_dot * dt + V_hat_ii(i-1);
end

%% Plots

% Plot de x

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t, V_hat, "LineWidth", 4)
plot(t, V_hat_ii, '-.',  "LineWidth", 2)
plot(t, Vm, '--', "LineWidth", 3)

legend("Estimado (i)", "Estimado (ii)", "Modelo de Referência")
xlabel("Tempo [s]")
ylabel("Velocidade")
title("Experiência 7(c) - Velocidade")
if saveFig
    saveas(gcf,'Exp7c_Velocidade.png')
end

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,e)
plot(t,e_ii)
legend("Estimado (i)", "Estimado (ii)")
xlabel("Tempo [s]")
ylabel("Erro")
title("Experiência 7(c) - Erro")
if saveFig
    saveas(gcf,'Exp7c_Erro.png')
end

% Plot de k

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,k)
plot(t,k_ii)
plot(t, k_star_vector, '--')
legend("Estimado (i)", "Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
ylim([-0.1, 1.5])
title("Experiência 7(c) - Parâmetro de controle k")
if saveFig
    saveas(gcf,'Exp7c_ParametroK.png')
end

% Plot de l

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l)
plot(t,l_ii)
plot(t, l_star_vector, '--')
legend("Estimado (i)","Estimado (ii)","Ideal")
xlabel("Tempo [s]")
title("Experiência 7(c) - Parâmetro de controle L")
if saveFig
    saveas(gcf,'Exp7c_ParametroL.png')
end

% Plot de delta

delta_star_vector = delta_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,delta, "LineWidth", 2)
plot(t,delta_ii, '--', "LineWidth", 2)
plot(t, delta_star_vector, '--')
legend("Estimado (i)","Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
title("Experiência 7(c) - Parâmetro de controle delta")
if saveFig
    saveas(gcf,'Exp7c_ParametroDelta.png')
end

%% Experiência 7 (c) - Indireto (i)

% Valores 
a = 0.5;
b = 1.5;
d = 10;

am = 0.5;
bm = 0.5;

% Valores ideais de controle
k_star = (am-a)/b;
l_star = bm/b;
delta_star = -d*a/b;

% Valores iniciais de adaptação
a0 = 1;
b0 = 1;
d0 = 1;
b_limite=1;
k0 = 1;
l0 = 1;
delta0 = 1;

gamma1 = 0.5;
gamma2 = 0.5;
gamma3 = 0.5;


% Valores de V para o controlador adaptativo
V_hat_indireto = zeros(1,length(t));
V0_hat = 0;
V_hat_indireto(1) = V0_hat;

% Vetores para progressão dos ganhos, do erro e da planta
k_indireto = zeros(1,length(t));
k_indireto(1) = k0;
l_indireto = zeros(1,length(t));
l_indireto(1) = l0;
delta_indireto = zeros(1,length(t));
delta_indireto(1) = delta0;
e_indireto = zeros(1,length(t));
e_indireto(1) = V_hat_indireto(1)-Vm(1);

a_indireto = zeros(1,length(t));
a_indireto(1) = a0;
b_indireto = zeros(1,length(t));
b_indireto(1) = b0;
h_indireto = zeros(1,length(t));
h_indireto(1) = d0;
u_indireto = zeros(1,length(t));

for i=2:length(t)
    % V com controlador adaptativo
    e_indireto(i-1) = V_hat_indireto(i-1) - Vm(i-1);
    
    % Adaptação dos parâmetros da planta
    a_indireto_dot = -gamma1*e_indireto(i-1)*V_hat_indireto(i-1);
    a_indireto(i) = a_indireto_dot*dt + a_indireto(i-1);
    
    if( abs(b_indireto(i-1)) > b_limite || ( (abs(b_indireto(i-1)) == b_limite) && (sign(e_indireto(i-1)*u_indireto(i-1)*b)) >= 0))
        b_indireto_dot = gamma2*e_indireto(i-1)*u_indireto(i);
    else
        b_indireto_dot = 0;
        if b_indireto(i-1) > 0
            b_indireto(i-1) = abs(b_limite);
        else
            b_indireto(i-1) = -abs(b_limite);
        end
    end
    b_indireto(i) = b_indireto_dot*dt + b_indireto(i-1);
    
    h_indireto_dot = -gamma3*e_indireto(i-1);
    h_indireto(i) = h_indireto_dot*dt + h_indireto(i-1);

    % Atualização da lei de controle
    
    k_indireto_i = (am-a_indireto(i))/b_indireto(i);
    k_indireto(i) = k_indireto_i;
    
    l_indireto_i = (bm)/b_indireto(i);
    l_indireto(i) = l_indireto_i;
    
    delta_indireto_i = -(h_indireto(i))/b_indireto(i);
    delta_indireto(i) = delta_indireto_i;
    
    u_hat_indireto = -k_indireto(i)*V_hat_indireto(i-1) + l_indireto(i)*Vs(i-1) + delta_indireto(i);
    u_indireto(i) = u_hat_indireto;

    % Simulação da planta
    
    V_hat_indireto_dot = -a*V_hat_indireto(i-1) + b*u_indireto(i) + d;
    V_hat_indireto(i) = V_hat_indireto_dot * dt + V_hat_indireto(i-1);
end

%% Experiência 7 (c) - Indireto (ii)

b = 1.5;

% Valores iniciais de adaptação
a0_ii = 1;
b0_ii = 1;
h0_ii = 1;
b_limite_ii=1;
k0_ii = 1;
l0_ii = 1;
delta0_ii = 1;

gamma1_ii = 0.5;
gamma2_ii = 0.5;
gamma3_ii = 0.5;


% Valores de V para o controlador adaptativo
V_hat_indireto_ii = zeros(1,length(t));
V0_hat_ii = 0;
V_hat_indireto_ii(1) = V0_hat_ii;

% Vetores para progressão dos ganhos, do erro e da planta
k_indireto_ii = zeros(1,length(t));
k_indireto_ii(1) = k0_ii;
l_indireto_ii = zeros(1,length(t));
l_indireto_ii(1) = l0_ii;
delta_indireto_ii = zeros(1,length(t));
delta_indireto_ii(1) = delta0_ii;
e_indireto_ii = zeros(1,length(t));
e_indireto_ii(1) = V_hat_indireto_ii(1)-Vm(1);

a_indireto_ii = zeros(1,length(t));
a_indireto_ii(1) = a0_ii;
b_indireto_ii = zeros(1,length(t));
b_indireto_ii(1) = b0_ii;
h_indireto_ii = zeros(1,length(t));
h_indireto_ii(1) = h0_ii;

a_planta = zeros(1,length(t));
a_planta(1) = 0.5 + 0.04/(1+V0_hat);
d_planta = zeros(1,length(t));
d_planta(1) = 0.2 + sin(0.02*0);


u_indireto_ii = zeros(1,length(t));

for i=2:length(t)
    % V com controlador adaptativo
    e_indireto_ii(i-1) = V_hat_indireto_ii(i-1) - Vm(i-1);
    
    % Adaptação dos parâmetros da planta
    a_indireto_dot = -gamma1_ii*e_indireto_ii(i-1)*V_hat_indireto_ii(i-1);
    a_indireto_ii(i) = a_indireto_dot*dt + a_indireto_ii(i-1);
    
    if( abs(b_indireto_ii(i-1)) > b_limite_ii || ( (abs(b_indireto_ii(i-1)) == b_limite_ii) && (sign(e_indireto_ii(i-1)*u_indireto_ii(i-1)*b)) >= 0))
        b_indireto_dot = gamma2_ii*e_indireto_ii(i-1)*u_indireto_ii(i);
    else
        b_indireto_dot = 0;
        if b_indireto_ii(i-1) > 0
            b_indireto_ii(i-1) = abs(b_limite_ii);
        else
            b_indireto_ii(i-1) = -abs(b_limite_ii);
        end
    end
    b_indireto_ii(i) = b_indireto_dot*dt + b_indireto_ii(i-1);
    
    h_indireto_dot = -gamma3_ii*e_indireto_ii(i-1);
    h_indireto_ii(i) = h_indireto_dot*dt + h_indireto_ii(i-1);

    % Atualização da lei de controle
    
    k_indireto_i2 = (am-a_indireto_ii(i))/b_indireto_ii(i);
    k_indireto_ii(i) = k_indireto_i2;
    
    l_indireto_i2 = (bm)/b_indireto_ii(i);
    l_indireto_ii(i) = l_indireto_i2;
    
    delta_indireto_i2 = -(h_indireto_ii(i))/b_indireto_ii(i);
    delta_indireto_ii(i) = delta_indireto_i2;
    
    u_hat_indireto_ii = -k_indireto_ii(i)*V_hat_indireto_ii(i-1) + l_indireto_ii(i)*Vs(i-1) + delta_indireto_ii(i);
    u_indireto_ii(i) = u_hat_indireto_ii;

    % Simulação da planta
    a_planta(i) = a + 0.04/(1+V_hat_ii(i-1));
    d_planta(i) = 0.2 + sin(0.02*t(i-1));
    
    V_hat_indireto_dot_ii = -a_planta(i)*V_hat_indireto_ii(i-1) + b*u_indireto_ii(i) + d_planta(i)*a_planta(i);
    V_hat_indireto_ii(i) = V_hat_indireto_dot_ii * dt + V_hat_indireto_ii(i-1);
end

%% Plot dos gráficos indireto

% Plot de velocidade

index_fig = index_fig + 1;
figure(index_fig)
hold on 
grid on
plot(V_hat_indireto)
plot(V_hat_indireto_ii)
plot(Vm)
legend("Estimado (i)", "Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
ylabel("Velocidade")
title("Velocidade pelo controlador indireto")
if saveFig
    saveas(gcf,'Exp7c_Velocidade_Indireto.png')
end

% Plot de e

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,e_indireto)
plot(t,e_indireto_ii)
legend("Estimado (i)", "Estimado (ii)")
xlabel("Tempo [s]")
ylabel("Erro")
title("Erro para os controladores indiretos")
if saveFig
    saveas(gcf,'Exp7c_Erro_Indireto.png')
end

% Plot de k

k_star_vector = k_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,k_indireto)
plot(t,k_indireto_ii)
plot(t, k_star_vector, '--')
legend("Estimado (i)", "Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
title("Parâmetro K do controlador indireto")
if saveFig
    saveas(gcf,'Exp7c_ParametroK_Indireto.png')
end

% Plot de l

l_star_vector = l_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,l_indireto)
plot(t,l_indireto_ii)
plot(t, l_star_vector, '--')
legend("Estimado (i)","Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
title("Parâmetro L do controlador indireto")
if saveFig
    saveas(gcf,'Exp7c_ParametroL_Indireto.png')
end

% Plot de delta

delta_star_vector = delta_star*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,delta_indireto)
plot(t,delta_indireto_ii)
plot(t, delta_star_vector, '--')
legend("Estimado (i)","Estimado (ii)", "Ideal")
xlabel("Tempo [s]")
title("Parâmetro \delta do controlador indireto")
if saveFig
    saveas(gcf,'Exp7c_ParametroDelta_Indireto.png')
end

% Plot de b

b_vector = b*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,b_indireto)
plot(t,b_indireto_ii)
plot(t, b_vector, '--')
legend("Estimado (i)","Estimado (ii)","Ideal")
xlabel("Tempo [s]")
title("Parâmetro b estimado da planta")
if saveFig
    saveas(gcf,'Exp7c_ParametroB_Indireto.png')
end

% Plot de a - i

a_vector = a*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,a_indireto)
plot(t, a_vector)
legend("Estimado","Ideal")
xlabel("Tempo [s]")
title("Parâmetro a estimado da planta - Caso (i)")
if saveFig
    saveas(gcf,'Exp7c_ParametroA_i_Indireto.png')
end

% Plot de a - ii

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,a_indireto_ii)
plot(t, a_planta)
legend("Estimado","Ideal")
xlabel("Tempo [s]")
title("Parâmetro a estimado da planta - Caso (ii)")
if saveFig
    saveas(gcf,'Exp7c_ParametroA_ii_Indireto.png')
end

% Plot de h - i

h_vector = (a*d)*ones(1, length(t));
index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,h_indireto)
plot(t, h_vector)
legend("Estimado","Ideal")
xlabel("Tempo [s]")
title("Parâmetro h estimado da planta - Caso (i)")
if saveFig
    saveas(gcf,'Exp7c_ParametroH_i_Indireto.png')
end

% Plot de h - ii

index_fig = index_fig + 1;
figure(index_fig)
hold on
grid on
plot(t,h_indireto_ii)
plot(t, d_planta.*a_planta)
legend("Estimado","Ideal")
xlabel("Tempo [s]")
title("Parâmetro h estimado da planta - Caso (ii)")
if saveFig
    saveas(gcf,'Exp7c_ParametroH_ii_Indireto.png')
end
