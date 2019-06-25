% Explore Sistina in vivo/ phantom cycling choices

%close all; 
clear 
N_cycles = 5; 
Amplitude_TQ = 0.3;
Amplitude_DQ = 0.5;
Amplitude_SQ = 1.0;
Q = [Amplitude_SQ, Amplitude_DQ, Amplitude_TQ];
tau1_s =10e-3;

%---------------------------------
%In vivo phase cycle
%---------------------------------
Xi1 = pi/2; 
Xi2 = -pi/2; 

%continuously sampled signal:
phi1_30 = pi/6:0.1:((2*pi)*N_cycles+pi/6);
phi1_210 = 7*(pi/6):0.1:((2*pi)*N_cycles+(pi/6)*7);

srate_cycle = 6;
ph_inc = (2*pi/srate_cycle);
phi1_6_30 = (pi/6):(ph_inc):(2*pi)*N_cycles+(pi/6)-1;
phi2_6_210 = 7*(pi/6):(ph_inc):((2*pi)*N_cycles+(pi/6)*7-1);



%---------------------------------
%Additional Phantom cycle
%---------------------------------
Xi0 = 0; 
%continuously
phi1_120 = 4*(pi/6):0.1:(2*pi)*N_cycles+4*(pi/6);
phi1_300 = 10*(pi/6):0.1:(2*pi)*N_cycles+10*(pi/6);
%6 step
phi3_6_120 = 4*(pi/6):(ph_inc):(2*pi)*N_cycles+4*(pi/6)-1;
phi4_6_300 = 10*(pi/6):(ph_inc):(2*pi)*N_cycles+10*(pi/6)-1;
%---------------------------------
%---------------------------------


for t = 1:1:50 
    TE_s(t) = (t-1)*1e-3; % 
    %TE_s(t) = 10e-3;
    
for o_Hz = 1:1:100 %B0 offset
    omega(o_Hz) =2*pi*(o_Hz-1);   
    
    %---------------------------------
    %In Vivo Phase cycle, Xi = ±pi/2:
    %---------------------------------
    % 'continuously' sampled signal
    [S1(:,o_Hz,t), S1_stim(:,o_Hz,t)]= Coherence_NasignalRelaxation(phi1_30,Xi1,omega(o_Hz),Q,TE_s(t),tau1_s); 
    [S2(:,o_Hz,t), S2_stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi1_210,Xi2,omega(o_Hz),Q,TE_s(t),tau1_s); 
    
    % 6 step sampled signal 
    [S1_6(:,o_Hz,t), S1_6stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi1_6_30,Xi1,omega(o_Hz),Q,TE_s(t),tau1_s);
    [S2_6(:,o_Hz,t), S2_6stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi2_6_210,Xi2,omega(o_Hz),Q,TE_s(t),tau1_s);
    
    
    %-----------------------------------
    %Additional Phantom cycle, Xi = 0 :
    %-----------------------------------
    % 'continuously' sampled signal 
    [S3(:,o_Hz,t), S3_stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi1_120,Xi0,omega(o_Hz),Q,TE_s(t),tau1_s); %start phase 120, Xi = pi/2
    [S4(:,o_Hz,t), S4_stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi1_300,Xi0,omega(o_Hz),Q,TE_s(t),tau1_s); %start phase 300, Xi = pi/2
    
    % 6 step sampled signal 
    [S3_6(:,o_Hz,t), S3_6stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi3_6_120,Xi0,omega(o_Hz),Q,TE_s(t),tau1_s);
    [S4_6(:,o_Hz,t), S4_6stim(:,o_Hz,t)] = Coherence_NasignalRelaxation(phi4_6_300,Xi0,omega(o_Hz),Q,TE_s(t),tau1_s);
end

end


% Spectra:
%-----------------------------------
Spec1 = fftshift(fft(S1),1)./size(S1,1); 
Spec2 = fftshift(fft(S2),1)./size(S2,1);

Spec3 = fftshift(fft(S3),1)./size(S3,1);
Spec4 = fftshift(fft(S4),1)./size(S4,1);

Spec1_6 = fftshift(fft(S1_6),1)./size(S1_6,1);
Spec2_6 = fftshift(fft(S2_6),1)./size(S2_6,1);

Spec3_6 = fftshift(fft(S3_6),1)./size(S3_6,1);
Spec4_6 = fftshift(fft(S4_6),1)./size(S4_6,1);


% Reco FFT
%-----------------------------------
Spec_invivo   = sqrt(Spec1_6.^2 + Spec2_6.^2);
Spec_phantom1  = sqrt(Spec1_6.^2 + Spec3_6.^2);
Spec_phantom2  = sqrt(Spec2_6.^2 + Spec4_6.^2);
Spec_all  = sqrt(Spec_phantom1.^2 + Spec_phantom2.^2);


% Alternative Reco from signal
%-----------------------------------
Signal_invivo  = sqrt((S1_6 + S2_6).^2); 
Signal_phantom1 = sqrt(abs(S1_6).^2 + abs(S3_6).^2); 
Signal_phantom2 = sqrt(abs(S2_6).^2 + abs(S4_6).^2);  
Signal_all  = sqrt(abs(Signal_phantom1).^2 + abs(Signal_phantom2).^2);

Spec_invivo_fromsignal = fftshift(fft((Signal_invivo)),1)./size(Signal_invivo,1); 
Spec_phantom_fromsignal = fftshift(fft((Signal_phantom1)),1)./size(Signal_phantom1,1);     
Spec_all_fromsignal = fftshift(fft((Signal_all)),1)./size(Signal_all,1); 


% compute frequency vector 
%for phase cycle:
nyq_cycle = srate_cycle/2;
pnts_cycle = size(S1_6,1);
f_cycle = linspace(-nyq_cycle,+nyq_cycle,floor(pnts_cycle));
%for ´"continuously" phase incremented data:
srate = size(phi1_30,2)/N_cycles;
nyq = srate/2;
pnts = size(S1,1);
f = linspace(-nyq,+nyq,floor(pnts));



%% Stim signal alone
Sstim_invivo    = sqrt((S1_stim + S2_stim).^2); %Sstim_invivo    = sqrt(abs(S1_stim).^2 + abs(S2_stim).^2);
Sstim_invivo_6  = sqrt((S1_6stim + S2_6stim).^2);

as(permute(S1_stim,[2 3 1]),'colormap','parula')
as(permute(S2_stim,[2 3 1]),'colormap','parula')
as(permute(S3_stim,[2 3 1]),'colormap','parula')
as(permute(S4_stim,[2 3 1]),'colormap','parula')

as(permute(Sstim_invivo_6,[2 3 1]),'colormap','parula')


%% show phasecycle for 2D data
%figure;hold on;plot(phi1_6_30,imag(S1_6stim),'.-');plot(phi1_6_30,imag(S2_6stim))
%figure;hold on;plot(phi1_30,imag(S1_stim),'.-');plot(phi1_30,imag(S2_stim))


%% Heatmaps 
as(permute(S1_6,[2 3 1]),'colormap','parula')
as(permute(S2_6,[2 3 1]),'colormap','parula')
as(permute(S3_6,[2 3 1]),'colormap','parula')
as(permute(S4_6,[2 3 1]),'colormap','parula')


%% Spectra
 
as(Spec1_6,'colormap','parula')
as(Spec2_6,'colormap','parula')

as(Spec_invivo,'colormap','parula')
as(Spec_phantom1,'colormap','parula')
as(Spec_phantom2,'colormap','parula')
as(Spec_all,'colormap','parula')








