clear all
close all
clc

% Generation
s = 2; %variance
nbp = 150; %N
A = 0.4; %magnitude A
Kmax = 10000;
PFA = zeros(1,Kmax); %probability of False alarm
PD = zeros(1,Kmax); %probability of detection
miteration = 1000; %iterations
PFA1 = 0.1; PFA2 = 0.01;

for m=1:miteration
    x0 = randn(1,nbp)*sqrt(s); %generation of hypothesis H0
    x1 = A + randn(1,nbp)*sqrt(s); %generation of hypothesis H1
    T0 = mean(x0); %computation of the statistics under H0
    T1 = mean(x1); %computation of the statistics under H1
    T00(m) = T0;
    T11(m) = T1;
    l = -0.5; %threshold
    step = 0.0015; %to optimize if necessary
    for k=1:Kmax
        seuil(k) = l;
        if T0>l PFA(k) = PFA(k)+1;end
        if T1>l PD(k) = PD(k)+1;end
        l = l+step;
    end
end
PFA = PFA/miteration;
PD = PD/miteration;
PFA_theorique = 1-normcdf(seuil,0,sqrt(s/nbp));
PD_theorique = 1-normcdf(seuil,A,sqrt(s/nbp));
figure
plot(seuil, PFA, 'Color', 'r')
hold on
plot(seuil, PFA_theorique, '--', 'Color', 'k')
plot(seuil, PD, 'Color', [0 0.5 0])
plot(seuil, PD_theorique, '--', 'Color', 'b')
legend('PFA Expérimentale','PFA Théorique', 'PD Expérimentale', 'PD Théorique')
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
histogram(T00,30,"Normalization","pdf")
hold on
xspan = -0.8:1e-4:0.5;
plot(xspan,normpdf(xspan,0,sqrt(s/nbp)))
legend("Histogramme", "Courbe théorique")
title("Statistique T(x) sous H_0")

figure
histogram(T11,30,"Normalization","pdf")
hold on
xspan = -0.2:1e-4:1;
plot(xspan,normpdf(xspan,A,sqrt(s/nbp)))
legend("Histogramme", "Courbe théorique")
title("Statistique T(x) sous H_1")

gamma1 = norminv(1-PFA1,0,sqrt(s/nbp));
gamma2 = norminv(1-PFA2,0,sqrt(s/nbp));

% Signal processing
load('signalconnu.mat','signalconnu');
x1 = signalconnu;
gamma = gamma2;
for n=1:length(x1)-nbp
    T(n) = mean(x1(n:n+nbp));
    if(T(n)>gamma)
        H(n) = 1;
    else
        H(n) = 0;
    end
end    
    
figure()
subplot(3,1,1);
plot(x1)
title('Signal avec A connu')
subplot(3,1,2);
plot(T)
hold on
yline(gamma)
title('Statistiques du signal avec A connu')
subplot(3,1,3);
plot(H, 'b', 'LineWidth', 2)
title('Détection avec PFA='+string(PFA2))
