%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Program to simulate the anomalous diffusion
%                         of phage in mucus via CTRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;                                                      %time step.
Sim_length = 1e+3;                           %no.of steps use 1e+3 for MSD.
J = 1e+3;                                %no. of trials. Use 10000 for MSD.

%system parameters
kBT=300*(1.38e-23);                                   %thermal temperature.
eta=.001;                                              %viscocity of water.
Rp = 60e-9;                                                  %phage radius.

gamma_0 = 6*pi*eta*Rp;                        %coeff. of friction in water.
D0 = kBT/gamma_0;

nu = linspace(1e-6,1.3,20); %linspace(1e-6,2,20);

A_data = [];
wait_min_data = [];

for kk=2:length(nu)

nu_k = nu(kk); %1e-50;
%t0 = dt*(gamma(1-nu_k))^(-1/nu_k);
t0 = dt*(sinc(nu_k))^(1/nu_k);                     %minimum diffusion time.

Mat1=zeros(J,2*Sim_length);                 %storage for particle location.

    for j =1:J
        x_old = 0;                                                  %x-i.c.                   
        y_old = 0;                                                  %y-i.c.
        Time =0;                                                 %time i.c.
        wt=0;                                            %waiting time i.c.
        
        %solve finite diff equation:
        for n=1:Sim_length

            if Time >= wt

                wx = rand - 0.5;
                wy = rand - 0.5;
               
                %Langevin x-dir.
                x_new = x_old  + sqrt(24*D0*dt)*wx;
                %Langevin y-dir.
                y_new = y_old  + sqrt(24*D0*dt)*wy;   
                
                x_old=x_new;
                y_old=y_new;
                
                %generate waiting time.  
                U = rand;
                tau_c = (t0)*(U^(-1/(nu_k)));
                wt = Time + tau_c - 0*dt;                      %waiting time.
                
           end

        Mat1(j,(2*n-1):(2*n))=[x_old y_old];
        Time = Time + dt;
        
        end
    end

%MSD computation of entire ensemble----------------------------------------
nData = Sim_length;                                    %no. of data points.
msd = zeros(nData-1,1);                                  %storage.
        for h=2:nData 
        particle_0 = 0*Mat1(:,1:2);                           %initial pts.  
        particle_h = Mat1(:,(2*h-1):2*h);                   %lag = h*tau.
        deltaCoords = particle_0 - particle_h;             %difference.
        squaredDisplacement = sum(deltaCoords.^2,2);            %dx^2+dy^2.
        msd(h-1) = mean(squaredDisplacement(:,1));    %MSD normalized.
        end
tau = 1:(nData-1);
tau = dt*tau';

%Linear Fit----------------------------------------------------------------
ii = 1:1:(nData-1);
tau1=tau(ii); msd1 = msd(ii);
x  = log(tau1);
y  = log(msd1);
p  = polyfit(x,y,1);

g = p(1);                                              %diffusion exponent.
lnD = p(2);                                    %log - diffusion coefficent.
Dfit = exp(lnD);
yfit = lnD + g*x;                              

%--------------------------------------------------------------------------
figure(kk) %plot msd w/ fitted msd
plot(log(tau),log(msd),'ko','MarkerSize',5,'LineWidth',2)
hold on
plot(log(tau),log(4*D0/(dt^(nu_k-1))) ...
     + nu_k*log(tau),'k--','LineWidth',2)
%plot(log(tau),log(Dfit) +  nu_k*log(tau),'k--','LineWidth',2)
xlabel('lag time (sec)')
ylabel('Ensemble-averaged MSD')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
hold off

figure(100)   %plot log-log msd w/ fitted line
plot(nu_k,g,'ko','MarkerSize',10,'LineWidth',2)
hold on
xlabel('Power law exp: \nu')
ylabel('Diffusion exp: \alpha')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

model = log(4*D0/((gamma(nu_k+1))*(dt^(nu_k-1)))) +  nu_k*log(tau);
err = log(msd)-model;

%{
figure(101)   %plot log-log msd w/ fitted line
hold on
plot(tau(1:1:end),kk + err(1:1:end),'ko','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k')
plot(tau(1:1:end),kk + 0*tau(1:1:end),'k','MarkerSize',10,'LineWidth',2)
xlabel('lag times')
ylabel('error')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow
%}

%minimum waiting time
wait_min = dt*(1/gamma(1-nu_k))^(1/nu_k);
A = Dfit;
A_data = [A_data; A];
wait_min_data = [ wait_min_data; wait_min];

figure(102)   %plot log-log msd w/ fitted line
hold on
plot(nu_k,A,'ko','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k')
xlabel('PowerLaw Exponent: \nu')
ylabel('Diffusion Coefficent: K')
axis square
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','bold')
drawnow

end
%write data to excel file
filename = 'MSD_data.xlsx';
C = [A_data wait_min_data nu'];
xlswrite(filename,C)
