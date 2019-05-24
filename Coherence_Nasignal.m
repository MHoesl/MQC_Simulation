function S=Coherence_Nasignal(phi1,xsi1,omega,Q,TE)

% parameters
%omega=2*pi*0;%rad/s imperfect B0
%phi1=30/180*pi;%rad
%phi2=phi1+pi/2;%rad
%TE=10e-3;%s
phi2=phi1+xsi1;
tau1=10e-3;%s
tau2=0.1e-3;%s

% amplitudes A of coherence pathways depend on (tau1, tau2, theta, t)
% dim 1 is p1=-1,+1
% dim 2 is p2=-3,-2,-1,+1,+2,+3
A=zeros(2,6);


% SQ pathways are set to 1/4
SQ=Q(1);     
A(1,3)=-SQ/4; %-1,-1
A(1,4)=-SQ/4; %-1,+1
A(2,3)=+SQ/4; %+1,-1
A(2,4)=+SQ/4; %+1,+1


% DQ pathawys are set to 1/4
DQ=Q(2);      
A(1,2)=-DQ/4; %-1,-2
A(1,5)=-DQ/4; %-1,+2
A(2,2)=+DQ/4; %+1,-2
A(2,5)=+DQ/4; %+1,+2


% TQ pathawys are set to 1/4
TQ=Q(3);      
A(1,1)=-TQ/4; %-1,-3
A(1,6)=-TQ/4; %-1,+3
A(2,1)=+TQ/4; %+1,-3
A(2,6)=+TQ/4; %+1,+3


S=myNasignal(phi1,phi2,omega,A,tau1,tau2,TE);

end

% Fleysher 2010 eq.[1] / eq.[11]
function s=myNasignal(phi1,phi2,omega,A,tau1,tau2,TE)

s=0;
p1=[-1 ,+1];
p2=[-3 ,-2 ,-1 ,+1 ,+2 ,+3];
for m=1:2
    for n=1:6
        s=s+exp(-1i*(p1(m)*phi1+(p2(n)-p1(m))*phi2)) .* B_na(p1(m),p2(n),tau1,tau2,omega,TE,A(m,n));
    end
end

end

% Fleysher 2010 eq.[2], A is the amplitude of the specified coherence pathway 
function B=B_na(p1,p2,tau1,tau2,omega,TE,A)

B = exp(-1i*(p1*tau1+p2*tau2)*omega) .* exp(1i*omega*TE) .* A;

end
        