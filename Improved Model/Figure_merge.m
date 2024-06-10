% Author: Yipeng Li (y9li@ucsd.edu)
% Load figures
fig1 = openfig('figure1b.fig', 'invisible');
ax1 = findall(fig1, 'type', 'axes');
fig2 = openfig('figure2b.fig', 'invisible');
ax2 = findall(fig2, 'type', 'axes');

% merged figure
figure;
n = length(ax1);

% colors
color_pcv_figure1 = [0 0 1]; % Solid blue for PCV from figure 1
color_pcv_figure2 = [0 0 1]; % Dashed blue for PCV from figure 2
color_vcv_figure1 = [1 0 0]; % Solid red for VCV from figure 1
color_vcv_figure2 = [1 0 0]; % Dashed red for VCV from figure 2

% plot
for i = 1:n
    subplot(2, 2, i);
    
    fig1_children = get(ax1(n - i + 1), 'Children'); 
    for j = 1:length(fig1_children)
        new_obj = copyobj(fig1_children(j), gca);
        if strcmp(new_obj.Type, 'line')
            if strcmp(new_obj.LineStyle, '--') 
                new_obj.Color = color_vcv_figure1; 
                new_obj.LineStyle = '-'; 
            else 
                new_obj.Color = color_pcv_figure1; % PCV solid
            end
        end
    end
 
    fig2_children = get(ax2(n - i + 1), 'Children');
    for j = 1:length(fig2_children)
        new_obj = copyobj(fig2_children(j), gca);
        if strcmp(new_obj.Type, 'line')
            if strcmp(new_obj.LineStyle, '--') % VCV dash
                new_obj.Color = color_vcv_figure2; % VCV dash
                new_obj.LineStyle = '--'; % dashed
            else % PCV solid
                new_obj.Color = color_pcv_figure2; % PCV dash 
                new_obj.LineStyle = '--'; % dashed
            end
        end
    end
    
    title(get(get(ax1(n - i + 1), 'Title'), 'String'));
    xlabel(get(get(ax1(n - i + 1), 'XLabel'), 'String'));
    ylabel(get(get(ax1(n - i + 1), 'YLabel'), 'String'));
    grid on;
    xlim([168 172]);
    xticks(168:172);
    xticklabels({'0', '1', '2', '3', '4'});
   
    if i == 2
        ylim([34 50]);
    elseif i == 3
        ylim([126 134]);
    elseif i == 4
        ylim([34 40]);
    end
end

subplot(2, 2, 2);
legend('PCV/CO_2(-)', 'VCV/CO_2(-)', 'PCV/CO_2(+)', 'VCV/CO_2(+)');
% subplot(2, 2, 2);
% legend('PCV Without CO2 Binding', 'VCV Without CO2 Binding', 'PCV With CO2 Binding', 'VCV With CO2 Binding');
% subplot(2, 2, 3);
% legend('PCV Without CO2 Binding', 'VCV Without CO2 Binding', 'PCV With CO2 Binding', 'VCV With CO2 Binding');
% subplot(2, 2, 4);
% legend('PCV Without CO2 Binding', 'VCV Without CO2 Binding', 'PCV With CO2 Binding', 'VCV With CO2 Binding');