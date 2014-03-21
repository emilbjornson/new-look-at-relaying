%This Matlab script can be used to generate Figure 8 of the paper:
%
%Emil Björnson, Michail Matthaiou, Mérouane Debbah, "A New Look at Dual-Hop
%Relaying: Performance Limits with Hardware Impairments" accepted for
%publication in IEEE Transcations on Communications.
%
%Download article: http://arxiv.org/pdf/1311.2634
%
%This is version 1.2 (Last edited: 2014-03-21).
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

clear all
close all


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));


%Relaying with transceiver hardware impairments


%Two scenarios: Rayleigh fading (1), Nakagami-m fading (2)
scenario = 2;

%Select fading parameters
if scenario == 1 %Rayleigh fading scenario
    alpha1 = 1; %Not used in Rayleigh fading distribution
    beta1 = 1;  %Variance at first hop
    
    alpha2 = 1; %Not used in Rayleigh fading distribution
    beta2 = 1;  %Variance at second hop
    
elseif scenario == 2 %Nakagami-m fading scenario
    
    alpha1 = 2; %alpha fading parameter at first hop
    beta1 = 1;  %beta fading parameter at first hop
    
    alpha2 = 2; %alpha fading parameter at second hop
    beta2 = 1;  %beta fading parameter at second hop
    
end


R = 4;   %Transmit information rate (per 2 channel uses)
x = 2.^R-1;  %Corresponding outage thresholds

N1 = 1;      %Normalized noise variance at relay
N2 = 1;      %Normalized noise variance at destination

SNR_dB = [20 30]; %SNR range at first hop
P1_dB = SNR_dB - 10*log10(alpha1*beta1); %Corresponding transmit power range at source (in dB)
P1 = 10.^(P1_dB/10); %Corresponding transmit power range at source

P2 = P1; %Set same transmit power range at relay.

%Range of Error Vector Magnitudes (EVMs) at first and second hop
kappaRange = 0:0.0001:0.3;


%Placeholders for analytic results: Amplify-and-Forward (AF) case
OP_ideal_af_f = zeros(length(kappaRange),length(SNR_dB)); %Outage probabilities, AF fixed gain, ideal hardware
OP_ideal_af_v = zeros(length(kappaRange),length(SNR_dB)); %Outage probabilities, AF variable gain, ideal hardware

OP_nonideal_af_f = zeros(length(kappaRange),length(SNR_dB)); %Outage probabilities, AF fixed gain, non-ideal hardware
OP_nonideal_af_v = zeros(length(kappaRange),length(SNR_dB)); %Outage probabilities, AF variable gain, non-ideal hardware

%Placeholders for analytic results: Decode-and-Forward (DF) case
OP_ideal_df = zeros(length(kappaRange),length(SNR_dB));    %Outage probabilities, DF, ideal hardware
OP_nonideal_df = zeros(length(kappaRange),length(SNR_dB)); %Outage probabilities, DF, non-ideal hardware


%Run simulations for different values of kappa1 and kappa2
for k = 1:length(kappaRange)
    
kappa1 = kappaRange(k); % Error Vector Magnitude of source
kappa2 = kappaRange(k); % Error Vector Magnitude of relay


    if scenario == 1
        
        Omega1 = beta1;    %New notation for channel variance at first hop
        Omega2 = beta2;    %New notation for channel variance at second hop
        
        
        %AF fixed gain relaying with ideal hardware. Outage probability is
        %given by [32]
        G_ideal_f = sqrt(P2./(P1*Omega1+N1));
        c_ideal_f = N2 ./ ( P1 .*(G_ideal_f.^2)*Omega1*Omega2 );
        
        OP_ideal_af_f(k,:) = 1 - 2 * exp(-N1*x./(P1*Omega1)) .* sqrt(c_ideal_f*x) .* besselk(1, 2*sqrt(c_ideal_f*x) );
        
        
        %AF variable gain relaying with ideal hardware. Outage probability is
        %given by [32]
        c_ideal_v = (N1*N2)./(P1.*P2*Omega1*Omega2);
        
        OP_ideal_af_v(k,:) = 1 - 2 * exp(-x*( N1./(P1*Omega1)+N2./(P2*Omega2))) .* sqrt(c_ideal_v*(x+x^2)) .* besselk(1, 2*sqrt(c_ideal_v*(x+x^2)) );
        
        
        %AF fixed gain relaying with non-ideal hardware. Outage probability is
        %given by Theorem 1
        G_nonideal_f = sqrt(P2./(P1*Omega1*(1+kappa1^2)+N1));
        A_nonideal_f = N2./(P1.* (G_nonideal_f.^2) *Omega1*Omega2);
        B_nonideal_f = N1*(1+kappa2^2)./(Omega1*P1);
        d = (kappa1^2+kappa2^2+kappa1^2*kappa2^2);
        
        if d*x>1
            OP_nonideal_af_f(k,:) = ones(size(B_nonideal_f));
        else
            OP_nonideal_af_f(k,:) = 1 - 2 * exp(-B_nonideal_f*x/(1-d*x)) .* sqrt(A_nonideal_f*x/(1-d*x)) .* besselk(1, 2*sqrt(A_nonideal_f*x/(1-d*x)) );
        end
        
        
        %AF variable gain relaying with non-ideal hardware. Outage probability
        %is given by Theorem 1
        A_nonideal_v = (N1*N2)./(P1.* P2 *Omega1*Omega2);
        B_nonideal_v = N1*(1+kappa2^2)./(Omega1*P1)+N2*(1+kappa1^2)./(Omega2*P2);
        
        if d*x>1
            OP_nonideal_af_v(k,:) = ones(size(B_nonideal_v));
        else
            OP_nonideal_af_v(k,:) = 1 - 2 * exp(-B_nonideal_v*x/(1-d*x)) .* sqrt(A_nonideal_v*(x+x^2))/(1-d*x) .* besselk(1, 2*sqrt(A_nonideal_v*(x+x^2))/(1-d*x) );
        end
        
        
        %DF relaying with ideal hardware. Outage probability is given in [33, Eq. (21)]
        prefactor1 = (N1*x/beta1);
        prefactor2 = (N2*x/beta2);
        OP_ideal_df(k,:) = 1 - exp(-prefactor1./P1 - prefactor2./P2);
        
        
        %DF relaying with non-ideal hardware. Outage probability is given by Theorem 2
        prefactor1 = (N1*x/(1-kappa1^2*x)/beta1);
        prefactor2 = (N2*x/(1-kappa2^2*x)/beta2);
        
        delta = max([kappa1^2 kappa2^2]);
        
        if delta*x>1
            OP_nonideal_df(k,:) = ones(size(B_nonideal_v));
        else
            OP_nonideal_df(k,:) = 1 - exp(-prefactor1./P1 - prefactor2./P2);
        end
        
    elseif scenario == 2
        
        %AF fixed gain relaying with ideal hardware. Outage probability is
        %given by [32]
        G_ideal_f = sqrt(P2./(P1*alpha1*beta1+N1));
        b1 = 0;
        b2 = N1./P1;
        c = N2./(P1.*G_ideal_f.^2);
        d = 0;
        
        OP_ideal_af_f(k,:) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF variable gain relaying with ideal hardware. Outage probability is
        %given by [32]
        b1 = N2./P2;
        b2 = N1./P1;
        c = N1*N2./(P1.*P2);
        d = 0;
        
        OP_ideal_af_v(k,:) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF fixed gain relaying with non-ideal hardware. Outage probability is
        %given by Theorem 1
        G_nonideal_f = sqrt(P2./(P1*alpha1*beta1*(1+kappa1^2)+N1));
        b1 = 0;
        b2 = (1+kappa2^2)*N1./P1;
        c = N2./(P1.*G_nonideal_f.^2);
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
        
        OP_nonideal_af_f(k,:) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %AF variable gain relaying with non-ideal hardware. Outage probability
        %is given by Theorem 1
        b1 = (1+kappa1^2)*N2./P2;
        b2 = (1+kappa2^2)*N1./P1;
        c = N1*N2./(P1.*P2);
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
        
        OP_nonideal_af_v(k,:) = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d);
        
        
        %DF relaying with ideal hardware. Outage probability is given in [33, Eq. (21)]
        prefactor1 = (N1*x/beta1);
        prefactor2 = (N2*x/beta2);
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        OP_ideal_df(k,:) = 1 - part1 .* part2;
        
        
        %DF relaying with non-ideal hardware. Outage probability is given by Theorem 2
        prefactor1 = (N1*x/(1-kappa1^2*x)/beta1);
        prefactor2 = (N2*x/(1-kappa2^2*x)/beta2);
        
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        delta = max([kappa1^2 kappa2^2]);
        
        if delta*x>1
            OP_nonideal_df(k,:) = ones(size(b1));
        else
            OP_nonideal_df(k,:) = 1 - part1 .* part2;
        end
        
    end

end


%Compute upper bounds on the OPs
upperKappaBoundAF = sqrt(sqrt(1/x+1)-1);
upperKappaBoundDF = 1/sqrt(x);


%Plot the results
figure; hold on; box on;

plot(kappaRange,OP_nonideal_af_f(:,1),'k','LineWidth',1);
plot(kappaRange,OP_nonideal_af_v(:,1),'b--','LineWidth',1);
plot(kappaRange,OP_nonideal_df(:,1),'r-.','LineWidth',1);

plot([upperKappaBoundAF upperKappaBoundAF],[1e-4 1],'k:','LineWidth',1);
plot([upperKappaBoundDF upperKappaBoundDF],[1e-4 1],'k:','LineWidth',1);

plot(kappaRange,OP_nonideal_af_f(:,2),'k','LineWidth',1);
plot(kappaRange,OP_nonideal_af_v(:,2),'b--','LineWidth',1);
plot(kappaRange,OP_nonideal_df(:,2),'r-.','LineWidth',1);

set(gca,'yscale','log')

legend('AF (Fixed Gain)','AF (Variable Gain)','DF','Necessary Condition','Location','SouthWest')

xlabel('Level of Impairments');
ylabel('Outage Probability (OP)');

axis([0 0.3 1e-4 1]);
