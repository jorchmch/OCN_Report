%% Proyecto: Theoretical and numerical dynamic of pulses in presence of
%% self phase modulation and negligible dispersion
%  SPM - SLOT 9, my Files is Class 5
%  SPM : SELF-PHASE MODULATION
%  GVD : GROUP VELOCITY DISPERSION

%%
% PROPAGATION METHOD FOR THE SOLUTION OF NLSE
% i Fz - beta2/2 Ftt + gamma*|F|^2 F=0
clc,clear, close all

% MATERIAL PROPERTIES
beta2=0; 
gamma=1;

% TEMPORAL COORDINATE
t0=-15;
t1=15;
nt=500;
t=linspace(t0,t1,nt);
deltat=(t1-t0)/(nt-1);

% SPATIAL COORDINATE
z0=0;
z1=10;
nz=500;
z=linspace(z0,z1,nz);
deltaz=(z1-z0)/(nz-1); % deltaz = h
h=deltaz;

% FREQUENCY COORDINATE 
indfreq=-nt/2:1:nt/2-1;
omega=(pi./t1).*indfreq;

% INPUT ENVELOPE

% pulseType =  'UNCHIRPED GAUSSIAN';
% tau0=1;
% amp=1;
% FIN=amp*exp(-t.^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE

% C=-2;
% amp= 1;
% tau0=1;
% FIN=amp*exp(-((1+1i*C)/2)*((t.^2)/(tau0^2))); %CHIRPED GAUSSIAN PULSE
% pulseType =  'CHIRPED GAUSSIAN';


% C=-2;
% tau0=1;
% FIN=sech(t/tau0).*exp(-1i*((C*t.^2)/(2*tau0^2))); % HYPERBOLIC CHIRPED SECH PULSES
% pulseType =  'HYPERBOLIC CHIRPED SECH';

% % SUPER GAUSSIAN PULSE
% m=3;
% C=5;
% tau0=1;
% FIN=exp(-(1+1i*C)/2*(t/tau0).^(2*m));
% pulseType =  'SUPER GAUSSIAN';


% 3 gaussian pulses
tau0=1;
amp=1;
pulseType =  '3 UGP';
tshift1=-10;
FIN1=amp*exp(-(t+tshift1).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 1
tshift2=0;
FIN2=1*amp*exp(-(t+tshift2).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 2
tshift3=10;
FIN3=amp*exp(-(t+tshift3).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 2
FIN=FIN1+FIN2+FIN3;

%FIN=amp*exp(-t.^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE


% FOURIER TRANSFORM - SIGNAL INPUT
FFIN=deltat*fftshift(fft(FIN));


%%%%%%%%%%%%%%%%%%%%%%
% CORE OF THE PROGRAM
%%%%%%%%%%%%%%%%%%%%%%

F=zeros(ceil(nt),floor(nz));
FF=zeros(ceil(nt),floor(nz));

F(:,1)=FIN;
FF(:,1)=FFIN;

q=FIN;

for loop_step=2:1:nz
   
    % LINEAR DISPERSIVE STEP
    qs=deltat*fftshift(fft(q)); %FFT
    qs_old=qs;
    
    prop=beta2/2*omega.^2;
    fact=1i*prop*h;
    qs=qs_old.*exp(fact); % calculation of the propagation in the frequency domain
    q=(1/deltat)*ifft(ifftshift(qs)); % coming back in the time domain
    
    %NONLINEAR CHI3 EFFECT STEP
    q_old=q;
    q=q_old.*exp(1i*gamma*abs(q_old).^2*h);
    
    %SAVE DATA EVERY NZ
    F(:,loop_step)=q;
    FF(:,loop_step)=deltat*fftshift(fft(q));
        
end


figure(1)
set(gcf,'Position',[100 100 900 300])
subplot(1,3,1)
mesh(t, z, abs(F)')
xlabel('Time (s)', 'FontSize', 12)
ylabel('Distance (z)', 'FontSize', 12)
zlabel('|F|', 'FontSize', 12)
title('Magnitude of F over Time and Space', 'FontSize', 12)
grid on


subplot(1,3,2)
mesh(omega, z, abs(FF)')
xlabel('Frequency (rad/s)', 'FontSize', 12)
ylabel('Distance (z)', 'FontSize', 12)
zlabel('|FFT F|', 'FontSize', 12)
title('Frequency Spectrum', 'FontSize', 12)
grid on
view(2)
xlim([-40 40])

subplot(1,3,3)
plot(omega, abs(FF(:,1)),'LineWidth',1.2)
hold on
plot(omega, abs(FF(:,end)),'LineWidth',1.2)
xlabel('Frequency (rad/s)', 'FontSize', 12)
ylabel('|FFT F|', 'FontSize', 12)
title('Frequency Spectrum:','FontSize', 12)
legend('Input', 'Output', 'FontSize', 10)
grid on
xlim([-40 40])

sgtitle(sprintf('%s Pulses',pulseType),'FontSize', 14)

% Figure 2: Phase Over Time and Space (Chirped Pulse)
figure(2)
set(gcf,'Position',[100 100 900 300])

subplot(1,3,1)
mesh(t, z, unwrap(angle(F))')
xlabel('Time (s)', 'FontSize', 12)
ylabel('Distance (z)', 'FontSize', 12)
zlabel('Phase (rad)', 'FontSize', 12)
title('Phase of F over Time and Space', 'FontSize', 12)
grid on
view(2)

subplot(1,3,2)
plot(t, unwrap(angle(F(:,1))),'LineWidth',1.2)
hold on
plot(t, unwrap(angle(F(:,end))),'LineWidth',1.2)
xlabel('Time (s)', 'FontSize', 12)
ylabel('Nonlinear Phase Shift \phi_{NL}', 'FontSize', 12)
title('Nonlinear Phase Shift', 'FontSize', 12)
legend('Input', 'Output', 'FontSize', 10)
grid on
% xlim([-5 5])

domegaIN=-gradient(unwrap(angle(F(:,1))),t);
domegaOUT=-gradient(unwrap(angle(F(:,end))),t);

subplot(1,3,3)
plot(t, domegaIN,'LineWidth',1.2)
hold on
plot(t, domegaOUT,'LineWidth',1.2)
xlabel('Time (s)', 'FontSize', 12)
ylabel('\delta \omega', 'FontSize', 12)
title('Induced Frequency Chirp', 'FontSize', 12)
legend('Input', 'Output', 'FontSize', 10)
grid on
% xlim([-5 5])





