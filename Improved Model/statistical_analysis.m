% Author: Yipeng Li (y9li@ucsd.edu)
close all
load('PCV_PA_NoCO2Binding.mat');
load('PCV_PA_WithCO2Binding.mat');
load('VCV_PA_NoCO2Binding.mat');
load('VCV_PA_WithCO2Binding.mat');

%load Po
PCV_Po = PCVno(:, 1);
PCV_Po_with = PCVwith(:, 1);
VCV_Po = VCVno(:, 1);
VCV_Po_with = VCVwith(:, 1);

% Performing t-test for Po
[~, p_value_PCV] = ttest2(PCV_Po, PCV_Po_with);
[~, p_value_VCV] = ttest2(VCV_Po, VCV_Po_with);
fprintf('Po literature vs. improved (PCV): P-value = %f\n', p_value_PCV);
fprintf('Po literature vs. improved (VCV): P-value = %f\n', p_value_VCV);

%load Pc
PCV_Pc = PCVno(:, 2);
PCV_Pc_with = PCVwith(:, 2);
VCV_Pc = VCVno(:, 2);
VCV_Pc_with = VCVwith(:, 2);
[~, p_value_PCVPc] = ttest2(PCV_Pc, PCV_Pc_with);
[~, p_value_VCVPc] = ttest2(VCV_Pc, VCV_Pc_with);
fprintf('Pc literature vs. improved (PCV): P-value = %f\n', p_value_PCVPc);
fprintf('Pc literature vs. improved (VCV): P-value = %f\n', p_value_VCVPc);

% %load Pao
% PCV_Pao = PCVno(:, 6);
% PCV_Pao_with = PCVwith(:, 6);
% VCV_Pao = VCVno(:, 6);
% VCV_Pao_with = VCVwith(:, 6);
% [~, PCV_Pao] = ttest2(PCV_Pao, PCV_Pao_with);
% [~, VCV_Pao] = ttest2(VCV_Pao, VCV_Pao_with);
% fprintf('P-value for Pao literature vs. improved (PCV): %f\n', PCV_Pao);
% fprintf('P-value for Pao literature vs. improved (VCV): %f\n', VCV_Pao);

% Perform t-test for Pao
[h_pcv_o2, p_pcv_o2] = ttest2(PCV_PA_O2_NoCO2Binding, PCV_PA_O2_withCO2Binding);
[h_vcv_o2, p_vcv_o2] = ttest2(VCV_PA_O2_NoCO2Binding, VCV_PA_O2_withCO2Binding);

% Perform t-test for Pac
[h_pcv_co2, p_pcv_co2] = ttest2(PCV_PA_CO2_NoCO2Binding, PCV_PA_CO2_withCO2Binding);
[h_vcv_co2, p_vcv_co2] = ttest2(VCV_PA_CO2_NoCO2Binding, VCV_PA_CO2_withCO2Binding);

% Display the t-test results with more precision
fprintf('Pao literature vs. improved (PCV): p-value = %.4e\n', p_pcv_o2);
fprintf('Pao literature vs. improved (VCV): p-value = %.4e\n', p_vcv_o2);
fprintf('Pac literature vs. improved (PCV): p-value = %.2e\n', p_pcv_co2);
fprintf('Pac literature vs. improved (VCV): p-value = %.2e\n', p_vcv_co2);

%% Calculate RMSE for Pac
% rmse_pcv_co2 = sqrt(mean((PCV_PA_CO2_NoCO2Binding - PCV_PA_CO2_withCO2Binding).^2));
% rmse_vcv_co2 = sqrt(mean((VCV_PA_CO2_NoCO2Binding - VCV_PA_CO2_withCO2Binding).^2));
% 
% % Calculate RMSE for Pao
% rmse_pcv_o2 = sqrt(mean((PCV_PA_O2_NoCO2Binding - PCV_PA_O2_withCO2Binding).^2));
% rmse_vcv_o2 = sqrt(mean((VCV_PA_O2_NoCO2Binding - VCV_PA_O2_withCO2Binding).^2));

% % Display the RMSE results
% fprintf('PCV CO2 RMSE: %f\n', rmse_pcv_co2);
% fprintf('VCV CO2 RMSE: %f\n', rmse_vcv_co2);
% fprintf('PCV O2 RMSE: %f\n', rmse_pcv_o2);
% fprintf('VCV O2 RMSE: %f\n', rmse_vcv_o2);

%% boxplot with significance star
alpha = 0.05;
% P-value from previous section
p_value_PCV = 0.589547;
p_value_VCV = 0.435807;
p_value_PCVPc = 0.000000;
p_value_VCVPc = 0.000000;
p_pcv_o2 = 0.060782;
p_vcv_o2 = 0.0067683;
p_pcv_co2 = 0.00;
p_vcv_co2 = 0.00;

% Function to add significance star
function addSignificanceStar(p, xpos, ypos)
    if p < 0.001
        text(xpos, ypos, '***', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    elseif p < 0.01
        text(xpos, ypos, '**', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    elseif p < 0.05
        text(xpos, ypos, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Plot for Po
figure;
subplot(2,2,1);
boxplot([PCV_Po, PCV_Po_with, VCV_Po, VCV_Po_with], 'Labels', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'});
set(gca, 'XTickLabel', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'}, 'TickLabelInterpreter', 'tex');
ylabel('P_o (mmHg)');
grid on;
text(0.01, 0.99, 'a', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
ymax = max([PCV_Po; PCV_Po_with; VCV_Po; VCV_Po_with]) * 1.1;
addSignificanceStar(p_value_PCV, 1.5, 120);
addSignificanceStar(p_value_VCV, 3.5, 120);

% Plot for Pc
subplot(2,2,2);
boxplot([PCV_Pc, PCV_Pc_with, VCV_Pc, VCV_Pc_with], 'Labels', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'});
set(gca, 'XTickLabel', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'}, 'TickLabelInterpreter', 'tex');
ylabel('P_c (mmHg)');
grid on;
text(0.01, 0.99, 'b', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
addSignificanceStar(p_value_PCVPc, 1.5, 36);
addSignificanceStar(p_value_VCVPc, 3.5, 36);

% Plot for Pao
subplot(2,2,3);
boxplot([PCV_PA_O2_NoCO2Binding, PCV_PA_O2_withCO2Binding, VCV_PA_O2_NoCO2Binding, VCV_PA_O2_withCO2Binding], 'Labels', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'});
set(gca, 'XTickLabel', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'}, 'TickLabelInterpreter', 'tex');
ylabel('P_a_o (mmHg)');
grid on;
text(0.01, 0.99, 'c', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
addSignificanceStar(p_pcv_o2, 1.5, 127);
addSignificanceStar(p_vcv_o2, 3.5, 127);

% Plot for Pac
subplot(2,2,4);
boxplot([PCV_PA_CO2_NoCO2Binding, PCV_PA_CO2_withCO2Binding, VCV_PA_CO2_NoCO2Binding, VCV_PA_CO2_withCO2Binding], 'Labels', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'});
set(gca, 'XTickLabel', {'PCV/CO_2(-)', 'PCV/CO_2(+)','VCV/CO_2(-)','VCV/CO_2(+)'}, 'TickLabelInterpreter', 'tex');
ylabel('P_a_c (mmHg)');
grid on;
text(0.01, 0.99, 'd', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
addSignificanceStar(p_pcv_co2, 1.5, 36.5);
addSignificanceStar(p_vcv_co2, 3.5, 36.5);



%% Error plot
% Calculate differences for Po, Pc, Pac, Pao
diff_pcv_po = PCV_Po_with - PCV_Po;
diff_vcv_po = VCV_Po_with - VCV_Po;
diff_pcv_pc = PCV_Pc_with - PCV_Pc;
diff_vcv_pc = VCV_Pc_with - VCV_Pc;
diff_pcv_co2 = PCV_PA_CO2_withCO2Binding - PCV_PA_CO2_NoCO2Binding;
diff_vcv_co2 = VCV_PA_CO2_withCO2Binding - VCV_PA_CO2_NoCO2Binding;
diff_pcv_o2 = PCV_PA_O2_withCO2Binding - PCV_PA_O2_NoCO2Binding;
diff_vcv_o2 = VCV_PA_O2_withCO2Binding - VCV_PA_O2_NoCO2Binding;

% Plot differences with RMSE line
figure;

% Po subplot
% subplot(2, 2, 1);
% hold on;
% plot(diff_pcv_po, 'b', 'DisplayName', 'PCV Po');
% plot(diff_vcv_po, 'r', 'DisplayName', 'VCV Po');
% yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero Line');
% xlabel('Sample Index');
% ylabel('Difference from P_o (mmHg)');
% grid on;
% ylim([-max(abs([diff_pcv_po; diff_vcv_po])) max(abs([diff_pcv_po; diff_vcv_po]))]);
% text(0.1, 0.9, 'a', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
% % Pc subplot
% subplot(2, 2, 2);
% hold on;
% plot(diff_pcv_pc, 'b', 'DisplayName', 'PCV Pc');
% plot(diff_vcv_pc, 'r', 'DisplayName', 'VCV Pc');
% yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero Line');
% xlabel('Sample Index');
% ylabel('Difference from P_c (mmHg)');
% grid on;
% ylim([-max(abs([diff_pcv_pc; diff_vcv_pc])) max(abs([diff_pcv_pc; diff_vcv_pc]))]);
% text(0.1, 0.9, 'b', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
% Pao subplot
subplot(1, 2, 1);
hold on;
plot(diff_pcv_o2, 'b', 'DisplayName', 'PCV O2');
plot(diff_vcv_o2, 'r', 'DisplayName', 'VCV O2');
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero Line');
xlabel('Sample Index');
ylabel('P_a_o Difference between Literature and Improved Model (mmHg)');
grid on;
ylim([-max(abs([diff_pcv_o2; diff_vcv_o2])) max(abs([diff_pcv_o2; diff_vcv_o2]))]);
text(0.01, 0.99, 'a', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
% Pac subplot
subplot(1, 2, 2);
hold on;
plot(diff_pcv_co2, 'b', 'DisplayName', 'PCV CO2');
plot(diff_vcv_co2, 'r', 'DisplayName', 'VCV CO2');
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero Line');
xlabel('Sample Index');
ylabel('P_a_c Difference between Literature and Improved Model (mmHg)');
legend('PCV (Improved)', 'VCV (Improved)','Literature Model (Zero Line)');
grid on;
ylim([-max(abs([diff_pcv_co2; diff_vcv_co2])) max(abs([diff_pcv_co2; diff_vcv_co2]))]);
hold off;
text(0.01, 0.99, 'b', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
%% Barplot
% mean_values = [
%     mean(PCV_PA_CO2_NoCO2Binding), mean(PCV_PA_CO2_withCO2Binding), mean(VCV_PA_CO2_NoCO2Binding), mean(VCV_PA_CO2_withCO2Binding);
%     mean(PCV_PA_O2_NoCO2Binding), mean(PCV_PA_O2_withCO2Binding), mean(VCV_PA_O2_NoCO2Binding), mean(VCV_PA_O2_withCO2Binding)
% ];
% 
% sem_values = [
%     std(PCV_PA_CO2_NoCO2Binding)/sqrt(length(PCV_PA_CO2_NoCO2Binding)), std(PCV_PA_CO2_withCO2Binding)/sqrt(length(PCV_PA_CO2_withCO2Binding)), std(VCV_PA_CO2_NoCO2Binding)/sqrt(length(VCV_PA_CO2_NoCO2Binding)), std(VCV_PA_CO2_withCO2Binding)/sqrt(length(VCV_PA_CO2_withCO2Binding));
%     std(PCV_PA_O2_NoCO2Binding)/sqrt(length(PCV_PA_O2_NoCO2Binding)), std(PCV_PA_O2_withCO2Binding)/sqrt(length(PCV_PA_O2_withCO2Binding)), std(VCV_PA_O2_NoCO2Binding)/sqrt(length(VCV_PA_O2_NoCO2Binding)), std(VCV_PA_O2_withCO2Binding)/sqrt(length(VCV_PA_O2_withCO2Binding))
% ];
% 
% % Bar plot with error bars
% figure;
% 
% % CO2 subplot
% subplot(1, 2, 1);
% bar(1:4, mean_values(1, :), 'FaceColor', 'flat');
% hold on;
% errorbar(1:4, mean_values(1, :), sem_values(1, :), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% title('CO2 Binding');
% xlabel('Condition');
% ylabel('PA CO2 (mmHg)');
% xticks(1:4);
% xticklabels({'PCV No CO2', 'PCV With CO2', 'VCV No CO2', 'VCV With CO2'});
% legend('Mean Values', 'Location', 'NorthEast');
% grid on;
% ylim([0 max(mean_values(1, :) + sem_values(1, :)) + 5]); % Adjust the y-axis limit for better visualization
% 
% % O2 subplot
% subplot(1, 2, 2);
% bar(1:4, mean_values(2, :), 'FaceColor', 'flat');
% hold on;
% errorbar(1:4, mean_values(2, :), sem_values(2, :), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% title('O2 Binding');
% xlabel('Condition');
% ylabel('PA O2 (mmHg)');
% xticks(1:4);
% xticklabels({'PCV No CO2', 'PCV With CO2', 'VCV No CO2', 'VCV With CO2'});
% legend('Mean Values', 'Location', 'NorthEast');
% grid on;
% ylim([0 max(mean_values(2, :) + sem_values(2, :)) + 5]); % Adjust the y-axis limit for better visualization
% 
% hold off;
% 
% % Prepare data for dot plot with error bars
% mean_values = [
%     mean(PCV_PA_CO2_NoCO2Binding), mean(PCV_PA_CO2_withCO2Binding), mean(VCV_PA_CO2_NoCO2Binding), mean(VCV_PA_CO2_withCO2Binding);
%     mean(PCV_PA_O2_NoCO2Binding), mean(PCV_PA_O2_withCO2Binding), mean(VCV_PA_O2_NoCO2Binding), mean(VCV_PA_O2_withCO2Binding)
% ];
% 
% sem_values = [
%     std(PCV_PA_CO2_NoCO2Binding)/sqrt(length(PCV_PA_CO2_NoCO2Binding)), std(PCV_PA_CO2_withCO2Binding)/sqrt(length(PCV_PA_CO2_withCO2Binding)), std(VCV_PA_CO2_NoCO2Binding)/sqrt(length(VCV_PA_CO2_NoCO2Binding)), std(VCV_PA_CO2_withCO2Binding)/sqrt(length(VCV_PA_CO2_withCO2Binding));
%     std(PCV_PA_O2_NoCO2Binding)/sqrt(length(PCV_PA_O2_NoCO2Binding)), std(PCV_PA_O2_withCO2Binding)/sqrt(length(PCV_PA_O2_withCO2Binding)), std(VCV_PA_O2_NoCO2Binding)/sqrt(length(VCV_PA_O2_NoCO2Binding)), std(VCV_PA_O2_withCO2Binding)/sqrt(length(VCV_PA_O2_withCO2Binding))
% ];
% 
%% Dot plot with error bars
% figure;
% 
% % CO2 subplot
% subplot(1, 2, 1);
% hold on;
% for i = 1:4
%     plot(i, mean_values(1, i), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
%     errorbar(i, mean_values(1, i), sem_values(1, i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% end
% title('CO2 Binding');
% xlabel('Condition');
% ylabel('PA CO2 (mmHg)');
% xticks(1:4);
% xticklabels({'PCV No CO2', 'PCV With CO2', 'VCV No CO2', 'VCV With CO2'});
% grid on;
% ylim([0 max(mean_values(1, :) + sem_values(1, :)) + 5]);
% 
% % O2 subplot
% subplot(1, 2, 2);
% hold on;
% for i = 1:4
%     plot(i, mean_values(2, i), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
%     errorbar(i, mean_values(2, i), sem_values(2, i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% end
% title('O2 Binding');
% xlabel('Condition');
% ylabel('PA O2 (mmHg)');
% xticks(1:4);
% xticklabels({'PCV No CO2', 'PCV With CO2', 'VCV No CO2', 'VCV With CO2'});
% grid on;
% ylim([0 max(mean_values(2, :) + sem_values(2, :)) + 5]);
% 
% hold off;