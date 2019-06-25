%% ----------------------------------------------------
%  simulation of the phase cycle for TQTPPI sequence
%  ----------------------------------------------------

%close all; 
%clear 

NCycles = 5; 
Amplitude_TQ = 0.3;
Amplitude_DQ = 0.5;
Amplitude_SQ = 1.0;
Q = [Amplitude_SQ, Amplitude_DQ, Amplitude_TQ];

tau1_s =10e-3;
ph0 = pi/2;
srate_cycle = 8;
ph_inc = (2*pi)/srate_cycle; 
phi1_cycle = ph0:ph_inc:((2*pi)*NCycles+ph_inc);

%'continuously' sampled signal
phi1 = ph0:0.1:((2*pi)* NCycles+ph_inc); 
srate = size(phi1,2)/NCycles;

Xiplus  = pi/2;
Ximinus = -pi/2; 


for o_Hz = 1:1:100 % B0 offset
    omega(o_Hz) =2*pi*(o_Hz-1); % offset in rad/s 

for t = 1:1:50 % increasing echo time
    TE_s(t) = (t-1)*1e-3; % 
    %TE(a) = 2e-3;

    S1(:,o_Hz,t) = Coherence_NasignalRelaxation(phi1,Xiplus,omega(o_Hz),Q,TE_s(t), tau1_s);
    S2(:,o_Hz,t) = Coherence_NasignalRelaxation(phi1,Ximinus,omega(o_Hz),Q,TE_s(t), tau1_s);
    
    S1_8(:,o_Hz,t) = Coherence_NasignalRelaxation(phi1_cycle,Xiplus,omega(o_Hz),Q,TE_s(t), tau1_s);
    S2_8(:,o_Hz,t)= Coherence_NasignalRelaxation(phi1_cycle,Ximinus,omega(o_Hz),Q,TE_s(t), tau1_s);

end
end


S_added = S1+ S2;
S8_added = S1_8+ S2_8;
 
% Spectra:
Spec1 = fftshift(fft(S1),1)/size(S1,1); 
Spec2 = fftshift(fft(S2),1)/size(S2,1);

Spec8_1 = fftshift(fft(S1_8),1)/size(S1_8,1);
Spec8_2 = fftshift(fft(S2_8),1)/size(S2_8,1);

Spec_added = fftshift(fft(S_added),1)/size(S_added,1); 
Spec8_added = fftshift(fft(S8_added),1)/size(S8_added,1);

%frequencies vector:
nyquist_bound_cycle = srate_cycle/2;
points_cycle = size(S8_added,1);
freq_vec = linspace(-nyquist_bound_cycle,+nyquist_bound_cycle,floor(points_cycle));


%frequencies vector "continuously"incremented data:
nyq = srate/2;
pnts = size(S_added,1);
f_vec = linspace(-nyq,+nyq,floor(pnts));



%% FOR 2D dataset:
% phase time figure
h1=figure(1); 
hold on
subplot(1,3,1); hold on; plot(phi1, imag(S1),'r');plot(phi1, imag(S2),'--r'); plot(phi1_cycle, imag(S1_8),'o','MarkerSize',3); plot(phi1_cycle, imag(S2_8),'ro','MarkerSize',3); title('Im \xi = ±\pi/2');
subplot(1,3,2); hold on; plot(phi1, real(S1),'b'); plot(phi1, real(S2),'--b');plot(phi1_cycle, real(S1_8),'o','MarkerSize',3);plot(phi1_cycle, real(S2_8),'bo','MarkerSize',3); title('Re \xi = ±\pi/2'); 
subplot(1,3,3); hold on; plot(phi1, imag(S_added),'r');plot(phi1, real(S_added),'--b','MarkerSize',3);plot(phi1_cycle, imag(S8_added),'ro','MarkerSize',3);plot(phi1_cycle, real(S8_added),'bo','MarkerSize',3);title('Added signals');
set(h1,'PaperSize',[12 6]);

% Spectra
h2=figure(2);
subplot(1,3,1);stem(abs(Spec1),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
subplot(1,3,2);stem(abs(Spec2),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
subplot(1,3,3);stem(abs(Spec_added),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
set(h2,'PaperSize',[12 6]);

h3=figure(3);
subplot(1,3,1);stem(abs(Spec8_1),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
subplot(1,3,2);stem(abs(Spec8_2),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
subplot(1,3,3);stem(abs(Spec8_added),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
set(h3,'PaperSize',[12 6]);


%% FOR 3D dataset:
% Heatmaps 
as(S1_8,'colormap','parula')
as(S2_8,'colormap','parula')
as(S8_added,'colormap','parula')

% Spectra
as(Spec8_1,'colormap','parula')
as(Spec8_2,'colormap','parula')
as(Spec8_added,'colormap','parula')

%frequency domain with useful frequency units showing nyquist limits
% h4=figure(4);
% subplot(1,3,1);stem3(abs(Spec1),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
% subplot(1,3,2);stem3(abs(Spec2),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
% subplot(1,3,3);stem3(abs(Spec_added),'ks-','linew',1,'markersize',2,'markerfacecolor','w')
% xlabel('frequency (1/cycle)'), ylabel('amplitude (a.u.)');
% set(h4,'PaperSize',[12 6]);



