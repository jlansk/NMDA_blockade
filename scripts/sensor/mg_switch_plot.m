%% Plots for understanding magnesium switch function
clear
close all

%% Select which figure you want (fig=1 or fig=2)
fig=2;

%% Plots

if fig == 1
    V=[-70:10:0]%fig 1;
    mg=linspace(-1, 1);% mg=[-1:0.1:1]
elseif fig == 2
    V=linspace(-70, 50);%fig2
    mg=[-1, -0.5, 0, 0.5, 1, 2, 4];% 2 3 4];%
end

alpha=exp(mg);
VN   =  10;                               % reversal Ca(NMDA)   

for v = 1:length(V)
    for a =1:length(alpha)        
        s(a, v) = 1.50265./(1 + 0.33*exp(-0.06.*alpha(a)'.*V(v)));% from mg_switch.m
    end
end


if fig == 1
    figure(1)
    % how s varies with mg
    for v = 1:length(V)
        plot(mg, s(:,v), 'LineWidth', 1)
        hold on
    end
    legend(string(V), 'location', 'NorthEastOutside')
    xlabel('blk_N_M_D_A')
    ylabel('m(V)')
    box off
    set(gcf, 'color', 'white')
end

if fig == 2
    % how m(V)/s varies with membrane potential (V)
    figure(2)
    for a = 1:length(mg)
        if mg(a) ~= 0
            plot(V, s(a,:), '-','LineWidth', 1)
            hold on
        elseif mg(a) ==0
            plot(V, s(a,:), '--', 'LineWidth', 1)
            hold on
        end
    end
    
    legend(string(mg), 'location', 'NorthEastOutside')
    xlabel('V'); xlim([-70 50])
    ylabel('m(V)')
    box off
    set(gcf, 'color', 'white')
    
end
