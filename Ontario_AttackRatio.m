% 0-11 years old
% 12+ years old

p =1-(0.03*0.57 + 0.86*0.80); %new vaccination coverage for Ontario as of Nov. 20, 2021

gamma1_0 = 1/7;
gamma2_0 = 1/7;

rho = 0.7;
theta = 0.0275;

% Calculation of baseline beta1_star for original strain
R = 2.5;
c11= 12.73; %Ontario

beta1_0 = (R*gamma1_0/c11); %this is the same as the baseline beta_1star in the manuscript

beta2_0 = 0.5*beta1_0; % kids half as susceptible

% Ontario
namecountry = 'Ontario';
N1 = 12932471;
N2 = 1801543;
S10=(1-(0.03*0.57 + 0.86*0.80))*N1; %susceptibles with the new vaccination coverage for Ontario as of Nov. 2, 2021
S20=N2;
E10=100;
E20=10;
I10=100;
I20=10;
c11 = 12.73; % (adult contact adult);
c12 = 1.03; % (adult contact kids);
c21 = 7.39; % (kids contact adults);
c22 = 4.3; % (kids contact kids);

red_contact = 1;

C = red_contact*[ c11, c12; c21, c22];

beta1_delta = beta1_0*1.4*1.6; %Nov. 2, 2021 for Ontario
beta2_delta = beta1_delta/2; 

% % Calculation of reproduction numbers
% Ra = beta1_delta/gamma1_delta_asy*c11*p;
% Rk = beta2_delta/gamma2_delta_asy*c22;
% 
% Re = 0.5*(Ra + Rk + sqrt( (Ra - Rk)^2 + 4*Ra*Rk*c12*c21 ) );

% Re = 0.5*(beta1*p*c11/gamma1 + beta2*c22/gamma2 + ...
%     sqrt( (beta1*p*c11/gamma1 - beta2*c22/gamma2)^2 + 4*beta1*p*c11/gamma1*beta2*c22/gamma2*c12*c21 ) );


% definition of baseline removal rates;
gamma1_delta = 1/7;
gamma2_delta = 1/7;

% gamma1_delta = 1/6;
% gamma2_delta = 1/6;

% gamma1_delta = 1/5;
% gamma2_delta = 1/5;

%gamma1_delta = 1/4;

% gamma1_delta = 1/3;
% gamma2_delta = 1/3;


%% Loop for final size as a function of gamma2, varying c22
gamma2_vector = 1./[0:0.1:10];
c22_vector = [c22, 0.5*c22, 0]; 

legendnames = {};

figure(1); hold on
for c22 = c22_vector
    
    C = red_contact*[ c11, c12; c21, c22];

    finalsize_kids_gamma2 = nan(length(gamma2_vector),1);

    for igamma = 1:length(gamma2_vector)
        
        gamma2_delta = gamma2_vector(igamma);

        Ra = beta1_delta/gamma1_delta*c11*p;
        Rk = beta2_delta/gamma2_delta*c22;
        Re = 0.5*(Ra + Rk + sqrt( (Ra - Rk)^2 + 4*Ra*Rk*c12*c21 ) );

        z = finalsize([beta1_delta;beta2_delta],[gamma1_delta;gamma2_delta],C,p,rho,E10,E20,N1,N2);

        if (Re > 1 && z(2)<=1 && z(2)>=0)
            finalsize_kids_gamma2(igamma) = 1-z(2);
        else
            finalsize_kids_gamma2(igamma) = nan;
            if (z(2)<0 || z(2)>1)
                display('Non admissible solution for final size')
            end
        end

    end
    
    legendnames{end+1} = ['$c_{22}=$',num2str(c22)];

    plot(1./gamma2_vector,finalsize_kids_gamma2,'LineWidth',3)
    
end

grid on
legend(legendnames,'Interpreter','latex','Location','NorthWest')
xlabel('$1/\gamma_2$ (days)', 'Interpreter','latex')
ylabel('$a_2$ (children attack ratio)', 'Interpreter','latex')
ylim([0 0.8]);
set(gca,'FontSize',25)



%% Loop for final size as a function of gamma2, varying beta2
gamma2_vector = 1./[0:0.1:10];
c22 = 4.3; % (kids contact kids) Ontario;
beta1_delta = beta1_0*1.4*1.6;
beta2_delta_vector = [0.7*beta1_delta, 0.6*beta1_delta, 0.5*beta1_delta, 0.4*beta1_delta]; %40-70% of beta1_delta 

legendnames = {};

figure(2); hold on
for ibeta2 = 1:length(beta2_delta_vector)
    
    C = red_contact*[ c11, c12; c21, c22];
    beta2_delta = beta2_delta_vector(ibeta2);

    finalsize_kids_gamma2beta = nan(length(gamma2_vector),1);

    for igamma2 = 1:length(gamma2_vector)
        
        gamma2_delta = gamma2_vector(igamma2);

        Ra = beta1_delta/gamma1_delta*c11*p;
        Rk = beta2_delta/gamma2_delta*c22;
        Re = 0.5*(Ra + Rk + sqrt( (Ra - Rk)^2 + 4*Ra*Rk*c12*c21 ) );

        z = finalsize([beta1_delta;beta2_delta],[gamma1_delta;gamma2_delta],C,p,rho,E10,E20,N1,N2);

        if (Re > 1 && z(2)<=1 && z(2)>=0)
            finalsize_kids_gamma2beta(igamma2) = 1-z(2);
        else
            finalsize_kids_gamma2beta(igamma2) = nan;
            if (z(2)<0 || z(2)>1)
                display('Non admissible solution for final size')
            end
        end

    end
    
    legendnames{end+1} = ['$\beta_{2}^{\star}=$',num2str(beta2_delta_vector(ibeta2),3)];
   

    plot(1./gamma2_vector,finalsize_kids_gamma2beta,'LineWidth',3)
    
end

grid on
legend(legendnames,'Interpreter','latex','Location','NorthWest')
xlabel('$1/\gamma_2$ (days)', 'Interpreter','latex')
ylabel('$a_2$ (children attack ratio)', 'Interpreter','latex')
ylim([0 0.8]);
set(gca,'FontSize',25)



%% Loop for final size as a function of c_{22}, varying gamma2
beta1_delta = beta1_0*1.4*1.6;
beta2_delta = beta1_delta/2;
gamma2_vector = 1./[3:7];
c22_vector = 0:0.1:6;

figure(3); hold on
legendnames = {};

for igamma = 1:length(gamma2_vector)

    gamma2_delta = gamma2_vector(igamma); 
    finalsize_kids_c22 = nan(length(c22_vector),1);

    for ic = 1:length(c22_vector)

        c22 = c22_vector(ic);
        C = red_contact*[ c11, c12; c21, c22];

        Ra = beta1_delta/gamma1_delta*c11*p;
        Rk = beta2_delta/gamma2_delta*c22;
        Re = 0.5*(Ra + Rk + sqrt( (Ra - Rk)^2 + 4*Ra*Rk*c12*c21 ) );

        z = finalsize([beta1_delta;beta2_delta],[gamma1_delta;gamma2_delta],C,p,rho,E10,E20,N1,N2);

        if (Re > 1 && z(2)<=1 && z(2)>=0)
            finalsize_kids_c22(ic) = 1-z(2);
        else
            finalsize_kids_c22(ic) = nan;
            if (z(2)<0 || z(2)>1)
                display('Non admissible solution for final size')
            end
        end

    end

    plot(c22_vector,finalsize_kids_c22,'LineWidth',3)
    
    legendnames{end+1} = ['$1/\gamma_2=$',num2str(1/gamma2_vector(igamma),1)];

end

grid on
legend(legendnames,'Interpreter','latex','Location','NorthWest') %NorthWest
xlabel('$c_{22}$', 'Interpreter','latex')
ylabel('$a_2$ (children attack ratio)', 'Interpreter','latex')
ylim([0 0.8]);
set(gca,'FontSize',25)



%% Loop for final size as a function of c_{22}, varying beta2
beta1_delta = beta1_0*1.4*1.6;
gamma2_delta = 1/7;

beta2_vector = [0.7*beta1_delta, 0.6*beta1_delta, 0.5*beta1_delta, 0.4*beta1_delta]; %40-70% of beta1_delta 
c22_vector = 0:0.1:6;

figure(4); hold on
legendnames = {};

for ibeta2 = 1:length(beta2_vector)

    beta2_delta = beta2_vector(ibeta2); 
    finalsize_kids_c22beta = nan(length(c22_vector),1);

    for ic = 1:length(c22_vector)

        c22 = c22_vector(ic);
        C = red_contact*[ c11, c12; c21, c22];

        Ra = beta1_delta/gamma1_delta*c11*p;
        Rk = beta2_delta/gamma2_delta*c22;
        Re = 0.5*(Ra + Rk + sqrt( (Ra - Rk)^2 + 4*Ra*Rk*c12*c21 ) );

        z = finalsize([beta1_delta;beta2_delta],[gamma1_delta;gamma2_delta],C,p,rho,E10,E20,N1,N2);

        if (Re > 1 && z(2)<=1 && z(2)>=0)
            finalsize_kids_c22beta(ic) = 1-z(2);
        else
            finalsize_kids_c22beta(ic) = nan;
            if (z(2)<0 || z(2)>1)
                display('Non admissible solution for final size')
            end
        end

    end

    plot(c22_vector,finalsize_kids_c22beta,'LineWidth',3)
    
    legendnames{end+1} = ['$\beta_2^{\star}=$',num2str(beta2_vector(ibeta2),3)];

end

grid on
legend(legendnames,'Interpreter','latex','Location','NorthWest') %NorthWest
xlabel('$c_{22}$', 'Interpreter','latex')
ylabel('$a_2$ (children attack ratio)', 'Interpreter','latex')
ylim([0 0.8]);
set(gca,'FontSize',25)



%% Definition of function to solve the final size equation 
function z = finalsize(betavector,gammavector,Cmatrix,coverage,rho,E10,E20,N1,N2)

    beta1 = betavector(1);
    beta2 = betavector(2);
    gamma1 = gammavector(1);
    gamma2 = gammavector(2);
   
    % without logarithm
    system = @(x) [ x(1) - (exp(beta1*rho*(Cmatrix(1,1)*coverage*x(1)/gamma1 + Cmatrix(1,2)*x(2)/gamma2) ) * exp(-beta1*rho*(Cmatrix(1,1)*coverage/gamma1 + Cmatrix(1,2)/gamma2))); ...
        x(2) - (exp(beta2*rho*(Cmatrix(2,1)*coverage*x(1)/gamma1 + Cmatrix(2,2)*x(2)/gamma2) ) * exp(-beta2*rho*(Cmatrix(2,1)*coverage/gamma1 + Cmatrix(2,2)/gamma2)))];

%     system = @(x) [ x(1) - (exp(beta1*rho*(Cmatrix(1,1)*coverage*x(1)/gamma1 + Cmatrix(1,2)*x(2)/gamma2) ) * exp(-beta1*rho*(Cmatrix(1,1)*coverage/gamma1 + Cmatrix(1,2)/gamma2)) * exp(-beta1*rho*(Cmatrix(1,1)*E10/(gamma1*N1) + Cmatrix(1,2)*E20/(gamma2*N2)))); ...
%         x(2) - (exp(beta2*rho*(Cmatrix(2,1)*coverage*x(1)/gamma1 + Cmatrix(2,2)*x(2)/gamma2) ) * exp(-beta2*rho*(Cmatrix(2,1)*coverage/gamma1 + Cmatrix(2,2)/gamma2)) * exp(-beta2*rho*(Cmatrix(2,1)*E10/(gamma1*N1) + Cmatrix(2,2)*E20/(gamma2*N2))))];

    init_vector = [0.1;0.4];
    z = fsolve(@(x) system(x), init_vector) % optimoptions('fsolve','Display','off')

end