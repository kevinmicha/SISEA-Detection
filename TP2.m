clear all
close all
clc
% Expressions are
% A_hat = mean(x)
% sigma0^2_hat = mean(x^2)
% sigma1^2_hat = mean((x-A_hat)^2)
% lambda = N*A^2/sigma^2
% E(x) = nu + lambda
% var(x) = 2nu + 4lambda

% Generation 
s = 2; %variance
nbp = 150; %N
A = 0.4; %magnitude A
lambda = nbp*A^2/s;
nu = 1;
Kmax = 1000;
PFA = zeros(1,Kmax); %probability of False alarm
PD = zeros(1,Kmax); %probability of detection
miteration = 1000; %iterations
PFA1 = 0.1; PFA2 = 0.01;

for m=1:miteration
    x0 = randn(1,nbp)*sqrt(s); %generation of hypothesis H0
    sigma0_squared_hat = mean(x0.^2);
    sigma1_squared_hat = mean((x0-mean(x0)).^2);
    T0 = nbp*log(sigma0_squared_hat/sigma1_squared_hat); %computation of the statistics under H0
    x1 = A + randn(1,nbp)*sqrt(s); %generation of hypothesis H1
    sigma0_squared_hat = mean(x1.^2);
    sigma1_squared_hat = mean((x1-mean(x1)).^2);
    T1 = nbp*log(sigma0_squared_hat/sigma1_squared_hat); %computation of the statistics under H1
    T00(m) = T0;
    T11(m) = T1;
    l = 0; %threshold
    step = 0.05; %to optimize if necessary
    for k=1:Kmax
        seuil(k) = l;
        if T0>l PFA(k) = PFA(k)+1;end
        if T1>l PD(k) = PD(k)+1;end
        l = l+step;
    end
end
PFA = PFA/miteration;
PFA_theorique = 1-chi2cdf(seuil,nu);
PD_theorique = 1-ncx2cdf(seuil,nu,lambda);
PD = PD/miteration;
figure
plot(seuil, PFA, 'Color', 'r')
hold on
plot(seuil, PFA_theorique, '--', 'Color', 'k')
plot(seuil, PD, 'Color', [0 0.5 0])
plot(seuil, PD_theorique, '--', 'Color', 'b')
legend('PFA Expérimentale', 'PFA Théorique', 'PD Expérimentale', 'PD Théorique')
xlabel('\gamma')
ylabel('Probabilité')
title('PFA(\gamma) et PD(\gamma)')

figure
axis([0,1.,0,1.])
hold on
plot(PFA,PD)
plot(PFA_theorique, PD_theorique)
grid on
xlabel('Probabilité de Fausse Alarme');
ylabel('Probabilité de Détection');
plot(PFA1,PD(find(PFA==PFA1+min(abs(PFA-PFA1)) | PFA==PFA1-min(abs(PFA-PFA1)),1)),'r*');
plot(PFA2,PD(find(PFA==PFA2+min(abs(PFA-PFA2)) | PFA==PFA2-min(abs(PFA-PFA2)),1)),'r*');
title('ROC')
legend('Expérimentale', 'Théorique')

figure
histogram(T00,50,"Normalization","pdf")
mean_h0 = mean(T00);
var_h0 = var(T00);
hold on
xspan = 0.05:1e-3:max(T00);
plot(xspan,chi2pdf(xspan,nu))
legend("Histogramme", "Courbe théorique")
title("Statistique T(x) sous H_0")

figure
histogram(T11,50,"Normalization","pdf")
mean_h1 = mean(T11);
var_h1 = var(T11);
hold on
xspan = 0.08:1e-3:max(T11);
plot(xspan,ncx2pdf(xspan,nu,lambda))
legend("Histogramme", "Courbe théorique")
title("Statistique T(x) sous H_1")

% Signal processing
sigma_1_rep = []; sigma_0_rep = []; A_rep = [];
load('Ainconnu.mat','x1');
gamma = chi2inv(1-PFA2,mean_h0);
for n=1:length(x1)-nbp
    A_hat = mean(x1(n:n+nbp));
    sigma0_squared_hat = mean(x1(n:n+nbp).^2);
    sigma1_squared_hat = mean((x1(n:n+nbp)-A_hat).^2);
    T(n) = nbp*log(sigma0_squared_hat/sigma1_squared_hat);
    if(T(n)>gamma)
        H(n) = 1;
        sigma_0_rep = [sigma_0_rep sigma0_squared_hat];
        sigma_1_rep = [sigma_1_rep sigma1_squared_hat];
        A_rep = [A_rep A_hat];
    else
        H(n) = 0;
    end
end    
    
figure()
subplot(3,1,1);
plot(x1)
title('Signal avec A inconnu')
subplot(3,1,2);
plot(T)
hold on
yline(gamma)
title('Statistiques du signal avec A inconnu')
subplot(3,1,3);
plot(H, 'b', 'LineWidth', 2)
title('Détection avec PFA='+string(PFA2))
