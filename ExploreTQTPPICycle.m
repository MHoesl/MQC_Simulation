%% 8 step TQTPPI phase cycling scheme 

clear 


NCycles = 1; 
Amplitude_TQ = 0.0;
Amplitude_DQ = 0.5;
Amplitude_SQ = 0.0;
Q = [Amplitude_SQ, Amplitude_DQ, Amplitude_TQ];

ph0 = pi/2;
srate_cycle = 8;
phase_inc = (2*pi)/srate_cycle; 
phi1_cycle = ph0:phase_inc:((2*pi)*NCycles+ph0-phase_inc);

%'continuously' sampled signal
phi1 = ph0:0.1:((2*pi)*NCycles+ph0); 
srate = size(phi1,2)/NCycles;


Xiplus = pi/2;
Ximinus = -pi/2; 

o_Hz = 1;
TE_s = 10e-3;
tau1_s = 10e-3;
for o_Hz = 1:2:50 %
   % TE = 2:2:30; % 
    omega(o_Hz) =2*pi*(o_Hz-1); 
    S1(:,o_Hz) = Coherence_Nasignal(phi1,Xiplus,omega(o_Hz),Q,TE_s);
    S2(:,o_Hz) = Coherence_Nasignal(phi1,Ximinus,omega(o_Hz),Q,TE_s);
    
    S1_8(:,o_Hz) = Coherence_Nasignal(phi1_cycle,Xiplus,omega(o_Hz),Q,TE_s);
    S2_8(:,o_Hz) = Coherence_Nasignal(phi1_cycle,Ximinus,omega(o_Hz),Q,TE_s);

end

S_added = S1(:,:) + S2(:,:);
S8_added = S1_8(:,:) + S2_8(:,:);


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
%freq_vec_halfspectrum = linspace(0,srate/2,floor(pnts/2)+1);


%frequencies vector "continuously"incremented data:
nyq = srate/2;
pnts = size(S_added,1);
f_vec = linspace(-nyq,+nyq,floor(pnts));





%% FIGURES

%phase time figure
h1=figure; 
hold on
subplot(1,3,1); hold on; plot(phi1, imag(S1),'r');plot(phi1, real(S1),'b'); plot(phi1_cycle, imag(S1_8),'ro'); plot(phi1_cycle, real(S1_8),'bo'); title('Im, Re, Xi = \pi/2');
subplot(1,3,2); hold on; plot(phi1, imag(S2),'r'); plot(phi1, real(S2),'b');plot(phi1_cycle, imag(S2_8),'ro');plot(phi1_cycle, real(S2_8),'bo'); title('Im, Re, Xi = -\pi/2'); 
subplot(1,3,3); hold on; plot(phi1, imag(S_added),'r');plot(phi1, real(S_added),'b');plot(phi1_cycle, imag(S8_added),'ro');plot(phi1_cycle, real(S8_added),'bo');title('Added signals Im, Re Xi = \pi/2');
set(h1,'PaperSize',[10 6]);



% frequency domain figure
h2=figure;
subplot(1,3,1); hold on; stem(freq_vec, abs(Spec8_1),'k','MarkerFaceColor',[1 0 0],'MarkerSize',3)  ;  title('Spec \xi = \pi/2');ylim([0 1]); xlim([-nyquist_bound_cycle-0.5 nyquist_bound_cycle+0.5]);xticks([-3 -2 -1 0 1 2 3 ]);
subplot(1,3,2); hold on; stem(freq_vec,abs(Spec8_2),'k','MarkerFaceColor',[1 0 0],'MarkerSize',3)    ; title('Spec \xi = -\pi/2');ylim([0 1]); xlim([-nyquist_bound_cycle-0.5 nyquist_bound_cycle+0.5]);xticks([-3 -2 -1 0 1 2 3 ]);
subplot(1,3,3); hold on; stem(freq_vec,abs(Spec8_added),'k','MarkerFaceColor',[1 0 0],'MarkerSize',3); title('Spec averaged signals'); ylim([0 1]); xlim([-nyquist_bound_cycle-0.5 nyquist_bound_cycle+0.5]);xticks([-3 -2 -1 0 1 2 3 ]);
set(h2,'PaperSize',[10 6]);



%% B0 inhomogeneity 
%phase time figure
[Phi1_8,Om] = meshgrid(phi1_cycle,omega/(2*pi));
[Phi,Om1] = meshgrid(phi1,omega/(2*pi));

h1=figure(1);
subplot(2,2,1); hold on; waterfall(Phi',Om1',imag(S1));stem3(Phi1_8',Om',imag(S1_8),'r'); title('Im Xi = +\pi/2'); xlim([0  12]);zlim([-2 2]); xticks([0 2 4 6 8 10 12]);view([40 -90 90]);
subplot(2,2,2); hold on; waterfall(Phi',Om1',real(S1));stem3(Phi1_8',Om',real(S1_8),'r'); title('Re Xi = +\pi/2'); xlim([0  12]);zlim([-0.002 0.002]); xticks([0 2 4 6 8 10 12]);view([40 -90 90]);
subplot(2,2,3); hold on; waterfall(Phi',Om1',imag(S2));stem3(Phi1_8',Om',real(S2_8),'r');title('Im Xi = -\pi/2');xlim([0  12]);zlim([-2 2]); xticks([0 2 4 6 8 10 12]);view([40 -90 90]);
subplot(2,2,4); hold on; waterfall(Phi',Om1',real(S2));stem3(Phi1_8',Om',real(S2_8),'r');title('Re Xi = -\pi/2'); xlim([0  12]);zlim([-0.002 0.002]); xticks([0 2 4 6 8 10 12]);view([40 -90 90]);
set(h1,'PaperSize',[8 7]);


%% Spectra
subplot(1,3,1); hold on; stem3(abs(Spec8_1)','MarkerSize',3);title('phase cycle Spec +');view([-30 45 45]);
subplot(1,3,2); hold on; stem3(abs(Spec8_2)','MarkerSize',3); title('phase cycle Spec -');view([-30 45 45]);
subplot(1,3,3); hold on; stem3(abs(Spec8_added)','MarkerSize',3); title('phase cycle av Spec');view([-30 45 45]);
set(h2,'PaperSize',[12 8]); %[width height]


