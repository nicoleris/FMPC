function [sim] = simulation(Tf, Ts, sys, mpc, x0)

T = ceil(Tf/Ts);
N = length(x0);

t = 0 : Ts : Tf;
x = zeros(N, T+1);
x(:, 1) = x0;

for i = 2:T+1
    u = mpc.get_u(x(:, i-1));
    x(:, i) = mpc.A * x(:, i-1) + mpc.B * u;
end

sim.t = t;
sim.x = x;

if N == 2
    if sys.StateName{1} == "vel_z"
        legend1 = "$\dot{z}$";
        data1 = "z velocity [m/s]";
    else
        legend1 = "$\dot{\gamma}$";
        data1 = "yaw velocity [rad/s]";
    end
    
    if sys.StateName{2} == "z"
        legend2 = "$z$";
        data2 = "z velocity [m]";
    else
        legend2 = "$\gamma$";
        data2 = "yaw [rad]";
    end
    
    figure;
    plot(t, x(1, :));
    hold on; grid on;
    plot(t, x(2, :));
    legend(legend1, legend2, 'interpreter','latex');
    xlabel('time [s]');
    ylabel(strcat(data1, " / ", data2), 'interpreter','latex');
    title('speed and position vs time');
      
elseif N == 4
    if sys.StateName{1} == "vel_pitch"
        legend1 = "$\dot{\alpha}$";
        data1 = "pitch velocity [rad/s]";
    else
        legend1 = "$\dot{\beta}$";
        data1 = "roll velocity [rad/s]";
    end
    
    if sys.StateName{2} == "pitch"
        legend2 = "$\alpha$";
        data2 = "pitch [rad]";
    else
        legend2 = "$\beta$";
        data2 = "roll [rad]";
    end
    
    if sys.StateName{3} == "vel_x"
        legend3 = "$\dot{x}$";
        data3 = "x velocity [m/s]";
    else
        legend3 = "$\dot{y}$";
        data3 = "y velocity [m/s]";
    end
    
    if sys.StateName{4} == "x"
        legend4 = "$x$";
        data4 = "x [m]";
    else
        legend4 = "$y$";
        data4 = "y [m]";
    end
    
    figure;

    subplot(2, 1, 1);
    plot(t, x(1, :));
    hold on; grid on;
    plot(t, x(2, :));
    legend(legend1, legend2, 'interpreter','latex');
    xlabel('time [s]');
    ylabel(strcat(data1, " / ", data2), 'interpreter','latex');
    title('Angular speed and position vs time');

    subplot(2, 1, 2);
    plot(t, x(3, :));
    hold on; grid on;
    plot(t, x(4, :));
    legend(legend3, legend4, 'interpreter','latex');
    xlabel('time [s]');
    ylabel(strcat(data3, " / ", data4), 'interpreter','latex');
    title('Axial speed and position vs time');

    sgtitle('Position and velocity vs time');
end
end


