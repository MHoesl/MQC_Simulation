function [S,S_stim,SQ_longitudinal]=Coherence_NasignalRelaxation(phi1,xsi1,omega,Q,TE,tau1)

% parameters

%turn ON/OFF stimulated echo
SQ_stim = 0;

phi2=phi1+xsi1;
tau2=0.1e-3;%s
%tau2=5e-3;%s
T2slow=30e-3;T2fast=4e-3;%s
exc_FA = pi/2;
T1=70e-3;%s

%Transverse Relaxation
SQ = Q(1) .*(3/5.*exp(-(TE+tau1)/T2slow) + 2/5.*exp(-(TE+tau1)/T2fast)).*sin(pi/2);
DQ = Q(2) .*(exp(-TE/T2slow) - exp(-TE/T2fast)).*(exp(-tau1/T2slow) - exp(-tau1/T2fast)).*exp(-tau2/T2slow).*sin(exc_FA).^5; 
TQ = Q(3) .*(exp(-TE/T2slow) - exp(-TE/T2fast)).*(exp(-tau1/T2slow) - exp(-tau1/T2fast)).*exp(-tau2/T2slow).*sin(exc_FA).^5;  


%Longitudinal Relaxation
SQ_l=1.0;
SQ_longitudinal = ((SQ_l.*(1-exp(-tau1/T1))).*(1-exp(-tau2/T1))).*(1-exp(-TE/T1));


% amplitudes A of coherence pathways depend on (tau1, tau2, theta, t)
% dim 1 is p1=-1,+1
% dim 2 is p2=-3,-2,-1,+1,+2,+3
A=zeros(2,6);


% SQ pathways are set to 1/4    
A(1,3)=-SQ/4; %-1,-1
A(1,4)=-SQ/4; %-1,+1
A(2,3)=+SQ/4; %+1,-1
A(2,4)=+SQ/4; %+1,+1

% DQ pathawys are set to 1/4    
A(1,2)=-DQ/4; %-1,-2
A(1,5)=-DQ/4; %-1,+3
A(2,2)=+DQ/4; %+1,-3
A(2,5)=+DQ/4; %+1,+3


% TQ pathawys are set to 1/4      
A(1,1)=-TQ/4; %-1,-3
A(1,6)=-TQ/4; %-1,+3
A(2,1)=+TQ/4; %+1,-3
A(2,6)=+TQ/4; %+1,+3


[S, S_stim]=myNasignal(phi1,phi2,omega,A,tau1,tau2,TE,SQ_stim,xsi1,T1,T2slow,T2fast);

end

% Fleysher 2010 eq.[1] / eq.[11]
function [s,s_stim]=myNasignal(phi1,phi2,omega,A,tau1,tau2,TE,SQ_stim,xsi1,T1,T2slow,T2fast)

s=0;
s_stim=0;
p1=[-1 ,+1];
p2=[-3 ,-2 ,-1 ,+1 ,+2 ,+3];
for m=1:2
    for n=1:6
        if ( abs(p2(n)) == 1 )
        s_stim = s_stim + B_stim(tau1,tau2,omega,TE,A(m,n),phi2,xsi1,T1,T2slow,T2fast).*SQ_stim;
        end
            
        s=s+ exp(-1i*(p1(m)*phi1+(p2(n)-p1(m))*phi2)) .* B_na(p1(m),p2(n),tau1,tau2,omega,TE,A(m,n));
        
    end  
end

s=s+s_stim;
end

% Fleysher 2010 eq.[2], A is the amplitude of the specified coherence pathway 
function B=B_na(p1,p2,tau1,tau2,omega,TE,A)

B = exp(-1i*(p1*tau1+p2*tau2)*omega) .* exp(1i*omega*TE) .* A;


end
        

% signal of stimulated echo at tau1 after the third pulse:
% add it to the transversal magnetization which is p2 = ±1
function B_stim = B_stim(tau1,tau2,omega,TE,A,phi2,xsi1,T1,T2slow,T2fast)
% approximation
T2starslow = T2slow; 
T2starfast = T2fast;


B_stim = (abs(A))./2 .* exp(-(tau2/T1)) .* cos(xsi1-omega*tau1).* exp(1i*omega*TE).*  (3/5.*exp(-(TE-tau1).^2./(2*T2starslow.^2)) + 2/5.*exp(-(TE-tau1).^2./(2*T2starfast.^2)));

end
