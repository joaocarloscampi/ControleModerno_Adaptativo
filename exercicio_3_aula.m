% Exercício 2 - Identificação de parâmetros

clc;
clear;
close all;

%% Interface 

disp("Entrada 1: Vetor inteiro");
disp("Entrada 2: Metade do vetor")

entrada = input("Qual você deseja? ");
clc;

if entrada==1
    half= false;
else 
    half= true;
end

disp("Entrada 1: Substituição direta");
disp("Entrada 2: Método do gradiente");
disp("Entrada 3: MMQ puro não recursivo");
disp("Entrada 4: MMQ puro recursivo");
disp("Entrada 5: MMQ não recursivo com fator de esquecimento");
disp("Entrada 6: MMQ recursivo com fator de esquecimento");

metodo = input("Qual você deseja? ");

%% Definição dos vetores

dt = 1;
t = 0:dt:14;
u = [1 0.8 0.6 0.4 0.2 0 0.2 0.4 0.6 0.8 1 0.8 0.6 0.4 0.2];
y = [0.9 2.5 2.4 1.3 1.2 0.8 0 0.9 1.4 1.9 2.3 2.4 2.3 1.3 1.2];

if half
    u = u(8:end);
    y = y(8:end);
    t = t(8:end);
end

%% a) Substituição direta

if metodo==1

theta = zeros(2,length(t)-1);

for i=2:length(t)-1
    M = [u(i), u(i-1);
         u(i+1), u(i)];
    y_i = [y(i); y(i+1)];
    
    theta(:, i) = M^(-1)*y_i;
end

end

%% b) Método do gradiente
if metodo==2 

% Vetores dos parâmetros
b0 = zeros(1, length(t)-1);
b1 = zeros(1, length(t)-1);
% Valores iniciais
b0(1) = 1;
b1(1) = 1;

theta = [b0; b1]';
phi = [u(2:end); u(1:(end-1))];

z = y(2:end);

gamma = 1;

z_hat = zeros(1, length(z));
e = zeros(1, length(z));

for i=1:length(z)-1
    z_hat(i) = theta(i,:)*phi(:,i);
    ms_2 = 1 + phi(:,i)'*phi(:,i);
    e(i) = (z(i) - z_hat(i))/ms_2;
    theta_dot = gamma*e(i)*phi(:,i)';
    theta(i+1,:) = theta_dot*dt + theta(i,:);
end

end
%% c) MMQ puro não recursivo
if metodo==3

% Vetores dos parâmetros
b0 = zeros(1, length(t)-1);
b1 = zeros(1, length(t)-1);
% Valores iniciais
b0(1) = 1;
b1(1) = 1;

%Parâmetros iniciais 
beta= 0;
Q0= [1 2; 2 1];
theta_0= [b0(1); b1(1)];

theta = zeros(2,length(t)-1);
phi = [u(2:end); u(1:(end-1))];

z = y(2:end);

%Inicialização
P= zeros(2,2,length(t)-1);
dt= 1;

for t_i = 1:length(t)-1
    integral= 0;
    
    for tal = 1:length(t)-1
        phi_i = phi(:,tal); % Valor de phi pro instante tal
        ms = 1+phi_i'*phi_i; % Valor de ms para o instante tal
        
        derivada = exp(-beta*(t_i-tal)) * (phi_i'*phi_i)/ms^(2);
        integral= integral + derivada*dt;
    end
    
    P(:,:,t_i) = (exp(-beta*t_i)*Q0 + integral)^(-1);
end    


for t_i = 1:length(t)-1
    integral_1= 0;
    
    for tal = 1:length(t)-1
        phi_i = phi(:,tal); % Valor de phi pro instante tal
        z_i = z(tal); % Valor de z pro instante tal
        ms = 1+phi_i'*phi_i; % Valor de ms para o instante tal
        
        derivada_1 = exp(-beta*(t_i-tal)) * (z_i*phi_i)/ms^(2);
        integral_1= integral_1 + derivada_1*dt;
    end
    
    theta(:,t_i) = P(:,:,t_i)*(exp(-beta*t_i)*Q0*theta_0 + integral);
end

end


%% d) MMQ puro recursivo

if metodo==4

% Vetores dos parâmetros
b0 = zeros(1, length(t)-1);
b1 = zeros(1, length(t)-1);
% Valores iniciais
b0(1) = 1;
b1(1) = 1;

theta = [b0; b1]';
phi = [u(2:end); u(1:(end-1))];

z = y(2:end);

%Valores iniciais de theta
theta_0= [b0(1); b1(1)];

% Matriz inicial P
P0 = [1 0 ;
      0 1];
  
% Progressão de theta
theta_v = [theta_0];
P_m= [];

for i=1:length(t)-1
    phi_i = phi(:,i); % Valor de phi pro instante i
    
    ms = 1+phi_i'*phi_i; % Valor de ms para o instante i
    
    P_d = -P0*phi_i*phi_i'*P0/(ms^2); % Calculo da derivada de P
    P = P0 + P_d*dt; % Atualização de P
    
    z_i = z(i);  % Valor de z para o instante i
    z_hat = theta_0'*phi_i; % Estimativa de z com os parametros atuais
    
    e = (z_i - z_hat)/ms^2; % Calculo de erro
    
    theta_d = P*e*phi_i;  % Calculo da derivada de theta
    theta = theta_0 + theta_d*dt; % Atualização de theta
    
    theta_v(:, end+1) = theta; % Armazenamento da progressao de theta_v
    
    % Valores para a proxima iteração
    P0 = P;
    theta_0 = theta;
end
theta= theta_v;
end

%% e) MMQ não recursivo com fator de esquecimento 
if metodo==5

% Vetores dos parâmetros
b0 = zeros(1, length(t)-1);
b1 = zeros(1, length(t)-1);
% Valores iniciais
b0(1) = 1;
b1(1) = 1;

%Parâmetros iniciais 
beta= 0.3;
Q0= [1 2; 2 1];
theta_0= [b0(1); b1(1)];

theta = zeros(2,length(t)-1);
phi = [u(2:end); u(1:(end-1))];

z = y(2:end);

%Inicialização
P= zeros(2,2,length(t)-1);
dt= 1;

for t_i = 1:length(t)-1
    integral= 0;
    
    for tal = 1:length(t)-1
        phi_i = phi(:,tal); % Valor de phi pro instante tal
        ms = 1+phi_i'*phi_i; % Valor de ms para o instante tal
        
        derivada = exp(-beta*(t_i-tal)) * (phi_i'*phi_i)/ms^(2);
        integral= integral + derivada*dt;
    end
    
    P(:,:,t_i) = (exp(-beta*t_i)*Q0 + integral)^(-1);
end    


for t_i = 1:length(t)-1
    integral_1= 0;
    
    for tal = 1:length(t)-1
        phi_i = phi(:,tal); % Valor de phi pro instante tal
        z_i = z(tal); % Valor de z pro instante tal
        ms = 1+phi_i'*phi_i; % Valor de ms para o instante tal
        
        derivada_1 = exp(-beta*(t_i-tal)) * (z_i*phi_i)/ms^(2);
        integral_1= integral_1 + derivada_1*dt;
    end
    
    theta(:,t_i) = P(:,:,t_i)*(exp(-beta*t_i)*Q0*theta_0 + integral);
end  

end

%% f) MMQ recursivo com fator de esquecimento

if metodo==6

% Vetores dos parâmetros
b0 = zeros(1, length(t)-1);
b1 = zeros(1, length(t)-1);
% Valores iniciais
b0(1) = 1;
b1(1) = 1;
beta= 0.3;

theta = [b0; b1]';
phi = [u(2:end); u(1:(end-1))];

z = y(2:end);

%Valores iniciais de theta
theta_0= [b0(1); b1(1)];

% Matriz inicial P
P0 = [1 0 ;
      0 1];
  
% Progressão de theta
theta_v = [theta_0];
P_m= [];

for i=1:length(t)-1
    phi_i = phi(:,i); % Valor de phi pro instante i
    
    ms = 1+phi_i'*phi_i; % Valor de ms para o instante i
    
    P_d = beta*P0 - P0*phi_i*phi_i'*P0/(ms^2); % Calculo da derivada de P
    P = P0 + P_d*dt; % Atualização de P
    
    z_i = z(i);  % Valor de z para o instante i
    z_hat = theta_0'*phi_i; % Estimativa de z com os parametros atuais
    
    e = (z_i - z_hat)/ms^2; % Calculo de erro
    
    theta_d = P*e*phi_i;  % Calculo da derivada de theta
    theta = theta_0 + theta_d*dt; % Atualização de theta
    
    theta_v(:, end+1) = theta; % Armazenamento da progressao de theta_v
    
    % Valores para a proxima iteração
    P0 = P;
    theta_0 = theta;
end
theta= theta_v;
end
disp("Os parâmetros encontrados são:")
disp(theta)