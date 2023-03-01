function plot_1d_profs(params, xlabels, titles, legendLabels, time) 
    %% plots the vertical 1d profiles in format (Z, params, data...)

    % Can set the screensize to keep figs consistent but doesn't work on
    % TACC
%     scrsz = get(0,'ScreenSize');
%     figure('Position', [0 0 scrsz(3)*0.9 scrsz(3)*0.9]);

    %figure()
    % figure for printing (keeps consistent sizing)
    figure('Renderer', 'painters', 'Position', [10 10 900 600])

    Z = params{size(params,2)};

    xlimits = zeros(8, 2);
    xlimits(1,:) = [-2, 2];
    xlimits(2,:) = [290, 360];
    xlimits(3,:) = [200,300];
    xlimits(4,:) = [0, 18];
    xlimits(5,:) = [0, 0.2];
    xlimits(6,:) = [0, 0.2];
    xlimits(7,:) = [0, 1];
    xlimits(8,:) = [0, 10^7];

    ylimits = [0,18];

    for i = 1:8
        subplot(2,4,i);
        for j = 1:size(params{i}, 1) % Loops over lines to plot
            P(i) = plot(params{i}(j,:), Z/10^3, 'LineWidth', 2); ...
                hold on; grid on;
        end
        legend(legendLabels, 'Location', 'best')
        xlabel(xlabels{i})
        ylabel('Height [km]')
        title(titles{i},'FontWeight','Normal')
        xlim([xlimits(i,1), xlimits(i,2)])
        ylim(ylimits)
    end
    sgtitle(strcat(num2str(time), ' hours'))
    exportgraphics(gcf,strcat('./plots/profiles/', 'profiles', '.gif'),...
        'Resolution',150, 'Append',true)

end
