% Example of FS-cell entering depolarization block (DB).
%  The idea: Increase reversal potential of potassium (EK) to enter DB.

clear
close all

%Set parameters.
T = 40000;              % # steps in simulation.
C = 1.0;                % Capacitance.
sigma = 0.5;            % Noise level.
gL  = 0.5;              % Leak conductance.
gNaF = 150;             % Fast sodium max conductance.
gKDR = 200;             % Fast potassium max conductance.
gCaH = 120;             % High-threhsold calcium max conductance.
gKM  = 320;             % M-current max conductance
gKv3 = 1200;            % Kv3.1 potassium max conductance.

%Define variables to save.
ic=0;                   % Set the initial conditions to 0.
V=[];                   % ... ***** ADD COMMENTS FOR EACH LINE ****
EK=[];                  % ...
I=[];
I0=ones(1,T)*-65;

%Start with EK at "baseline" level.
EK0=ones(1,T)*-100;
[V0,t,mNaF,hNaF,mKDR,mCaH,kV,mKM,ic] = traub_edit(T, I0, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic);
V = [V;V0];
I = [I,I0];
EK= [EK,EK0];

%Then, increase EK.
EK0=ones(1,T)*-60;
[V0,t,mNaF,hNaF,mKDR,mCaH,kV,mKM,ic] = traub_edit(T, I0, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic);
V = [V;V0];
I = [I,I0];
EK= [EK,EK0];
Ki = 130;
Ko = Ki*exp(EK/26.64);
t = (1:length(V))*(t(10)-t(9));

%Plot the results.
figure(10)
clf
set(gcf, 'Position', [0, 500, 500, 150])
set(gca,'FontSize', 12)
plot(t,V, 'k', 'LineWidth',2)
axis tight
ylim([-100 50])
set(gca, 'box','off')
ylabel('Voltage [mV]')
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'XTick', [], ...
           'Color','none',...
           'XColor','k','YColor','k', ...
           'FontSize', 12);
hl2 = line(t,Ko,'Color','r','LineWidth',2,'Parent',ax2);
axis tight
ylim([2 15])
ylabel('K_o [mM]')

