%% FELYSHER and RELAXATION:


close all; 
clear 

Amplitude_TQ = 0.3;
Amplitude_DQ = 0.5;
Amplitude_SQ = 1.0;
Q = [Amplitude_SQ, Amplitude_DQ, Amplitude_TQ];

N_cycles = 5; 
tau1 = 10e-3;
Xipi2 = -pi/2;   
Xi0 = 0;

phi10 = 0:0.1:(2*pi)*N_cycles;  
phi1 =pi/2:0.1:((2*pi))*N_cycles+pi/2;


samplingrate = size(phi1,2)/N_cycles;
samplingrate_cycle = 6; 
ph_inc = (2*pi)/samplingrate_cycle; 

phi_Xipi2 = (pi/2):(ph_inc):((2*pi))*N_cycles+pi/2 -1; % belongs to Xi = pi/2, the original phase cycle
phi_Xi0 = (0):(ph_inc):((2*pi))*N_cycles-1; % belongs to Xi = 0, the additional phase cycle


omega = 2*pi*0;% offset in rad/s 

for o_Hz = 1:1:100 % B0 offset
    omega(o_Hz) =2*pi*(o_Hz-1); % offset in rad/s 
    %omega(o_Hz) =2*pi*25; %(o-1); 


for t = 1:1:50
    TE_s(t) = (t-1) * 1e-3; % 
    %TE(a) = 10e-3;

    [S_Xipi2(:,o_Hz,t),S_stimXipi2(:,o_Hz,t)]   = Coherence_NasignalRelaxation(phi1,Xipi2,omega(o_Hz),Q,TE_s(t), tau1);
    [S_Xi0(:,o_Hz,t) ,S_stimXi0(:,o_Hz,t)]      = Coherence_NasignalRelaxation(phi10,Xi0,omega(o_Hz),Q,TE_s(t), tau1);
    [S_Xipi2_6(:,o_Hz,t),S_stimXipi2_6(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi_Xipi2,Xipi2,omega(o_Hz),Q,TE_s(t), tau1);
    [S_Xi0_6(:,o_Hz,t),S_stimXi0_6(:,o_Hz,t)]   = Coherence_NasignalRelaxation(phi_Xi0,Xi0,omega(o_Hz),Q,TE_s(t), tau1);
    
end
end

% Spectra
SpecXipi2 = fftshift(fft(S_Xipi2),1)/size(S_Xipi2,1);
SpecXi0 = fftshift(fft(S_Xi0),1)/size(S_Xi0,1); 
SpecXipi2_6 = fftshift(fft(S_Xipi2_6),1)/size(S_Xipi2_6,1);
SpecXi0_6 = fftshift(fft(S_Xi0_6),1)/size(S_Xi0_6,1);


% Added spectra:
Spec_plus   = 1/2*(SpecXi0 + 1i*SpecXipi2);
Spec_minus  = 1/2*(SpecXi0 - 1i*SpecXipi2);
Spec12_plus = 1/2*(SpecXi0_6 + 1i*SpecXipi2_6);
Spec12_minus= 1/2*(SpecXi0_6 - 1i*SpecXipi2_6); 
%sum of squares reco
Spec12 = sqrt(Spec12_plus.^2 + Spec12_minus.^2);


%sum of squares reco
Spec_sos_reco  = sqrt(SpecXi0_6.^2 + SpecXipi2_6.^2);


%% Recover signal with known B0 offset:
%tau1 = 10e-3;
%Spectrum = (S_plus.*exp(+1i*omega*tau1) - S_minus.*exp(-1i*omega*tau1)) .*exp(-1i*omega*TE);
%Spectrum_12step = (S12_plus.*exp(+1i*omega*tau1) - S12_minus.*exp(-1i*omega*tau1)) .*exp(-1i*omega*TE);

%% show Signal
as(S_Xipi2,'ColorMap',parula)
as(S_Xipi0,'ColorMap',parula)


%% show Spectra before reco

as(SpecXipi2_6,'ColorMap',parula)
as(SpecXi0_6,'ColorMap',parula)
as(Spec12_plus,'ColorMap',parula)
as(Spec12_minus,'ColorMap',parula)


%% Recover signal with known B0 offset:
tau1 = 10e-3;
%a = 11;
for t = 1:1:50
    TE_s(t) = (t-1) * 1e-3; % 
for o_Hz = 1:1:100
    omega(o_Hz) =2*pi*(o_Hz-1); 

    Spectrum(:,o_Hz,t) = (Spec_plus(:,o_Hz,t).*exp(+1i*omega(o_Hz)*tau1) - Spec_minus(:,o_Hz,t).*exp(-1i*omega(o_Hz)*tau1)) .*exp(-1i*omega(o_Hz)*TE_s(t));
    Spectrum_12step(:,o_Hz,t) = (Spec12_plus(:,o_Hz,t).*exp(+1i*omega(o_Hz)*tau1) - Spec12_minus(:,o_Hz,t).*exp(-1i*omega(o_Hz)*tau1)) .*exp(-1i*omega(o_Hz)*TE_s(t));
end
end

as(Spectrum_12step,'ColorMap',parula)

%compare to 
as(Spec12,'ColorMap',parula)


%% compute frequencies vector:
nyq_cycle = samplingrate_cycle/2; pnts_cycle = size(S_Xi0_6,1);
f_cycle = linspace(-nyq_cycle,+nyq_cycle,floor(pnts_cycle));
%for showing only half of the spectrum f = linspace(0,srate/2,floor(pnts/2)+1);
nyq = samplingrate/2; pnts = size(S_Xi0,1); f = linspace(-nyq,+nyq,floor(pnts));


