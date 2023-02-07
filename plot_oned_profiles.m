function plot_oned_profiles(Z, U,TH, QT, QC, TKE, TKE_HOR, TKE_W, WAT_FLUX, Z_CP, U_CP,TH_CP, QT_CP, QC_CP, TKE_CP, TKE_HOR_CP, TKE_W_CP, WAT_FLUX_CP)

%% Make the plots 
figure()
subplot(241);
plot(U_CP, Z/10^3,'LineWidth',2);hold on
plot(U, Z/10^3,'LineWidth',2);hold on;grid on
legend('CP', 'NO CP')
xlabel("u (ms^{-1})")
ylabel('Height [km]')

subplot(242);
plot(TH_CP,Z/10^3,'LineWidth',2); hold on
plot(TH,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("\theta (K)")

subplot(243);
plot(QT_CP,Z/10^3,'LineWidth',2); hold on
plot(QT,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("q_t (gkg^{-1})")

subplot(244);
plot(QC_CP,Z/10^3,'LineWidth',2);hold on
plot(QC,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("q_c (gkg^{-1})")

subplot(245);
plot(TKE_CP,Z/10^3,'LineWidth',2); hold on
plot(TKE,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("TKE (m^2s^{-2})")

subplot(246);
plot(TKE_HOR_CP,Z/10^3,'LineWidth',2);hold on
plot(TKE_HOR,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("1/2(u'^2+v'^2) (m^2s^{-2})")

subplot(247);
plot(TKE_W_CP,Z/10^3,'LineWidth',2); hold on
plot(TKE_W,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("1/2(w'^2) (m^2s^{-2})")

subplot(248);
plot(WAT_FLUX_CP.*10^2,Z/10^3,'LineWidth',2); hold on
plot(WAT_FLUX.*10^2,Z/10^3,'LineWidth',2);grid on;
legend('CP', 'NO CP')
ylabel('Height [km]')
xlabel("\times 10^{-2} (w'q_t') (mgs^{-1}kg^{-1})")