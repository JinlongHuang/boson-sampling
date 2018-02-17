% Plot for time efficiency
% ©2018 Jin-Long Huang All Right Reserved
% -------------------------------------------------------------------------

x = [6,7,8,9,10];
y1= [log10(7),log10(17),log10(2748),log10(74160),log10(2393280)];
y2= [log10(0.7),log10(73),log10(92),log10(1011),log10(1165)];
plot(x,y1, '-*', x, y2, '-^','LineWidth', 1.5)
ax = gca;
ax.FontSize = 18;
title('Time Efficiency');
xlabel('# Photons','FontSize',18);
ylabel('Log(Time/s)','FontSize',18);
legend({'Glynn','Gurvits'},'FontSize',16);