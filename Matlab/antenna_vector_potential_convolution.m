profile on
tic
%% parameters
load("parameters.mat","f","dr")
%dr = '..';

cpw_pos = [250:1:259 272:1:281 295:1:304]; % position of CPW lines

T = 200; % number of timesteps to use
Ms = 1.59E+5;
mu0 = 4*pi*1e-7;

metal = 160E-9; % metallization thickness

%% definitions
n_cpw = length(cpw_pos);

V_real(1:length(f)) = 0;
V_imag(1:length(f)) = 0;

[Mx{T}, My{T}, Mz{T}, Mx_pad, My_pad, Mz_pad]  = deal([]); 

for i_f=1:length(f) % loop through frequencies

    disp(['Processing #',num2str(i_f),'/',num2str(length(f)),' (f=',num2str(f(i_f)),' GHz)'])

    load([dr,'/results_compressed_f',num2str(i_f)],"Mx","My","Mz","dx","dy","dz","I","J","K","time")
    time = time(1:T);

    A_ = Ms*mu0/(4*pi)*dx*dy*dz;    % multiplicative factor in vector potential

    %% create coordinates for curl
    [X,Y,Z] = meshgrid(((1:I)*dx),((1:J)*dy),((1:K+1)*dz));

    %% calculate distance kernel in the plane above YIG, restricted to CPW
    x_ = (1:I+cpw_pos(end)-cpw_pos(1)+1)*dx;
    y_ = (1:2*J)*dy;
    dx2 = (x_-x_(I-cpw_pos(1)+1)).^2;
    dy2 = (y_-y_(J)).^2;
    [DX2,DY2] = meshgrid(dx2,dy2);

    distance_kernel_c = 1./sqrt(DX2 + DY2 + (dz/2 + metal/2).^2);   % center
    distance_kernel_t = 1./sqrt(DX2 + DY2 + (metal/2).^2);          % top surface
    distance_kernel_b = 1./sqrt(DX2 + DY2 + (dz + metal/2).^2);     % bottom surface

    nbytes = fprintf('Timestep 0 of %d', T);
    for i_t = 1:T   % loop thourgh timesteps
        %% print progress
        fprintf(repmat('\b',1,nbytes))
        nbytes = fprintf('Timestep %d of %d\n', i_t, T);

        %% pad for curl
        Mx_pad = (padarray(Mx{i_t},[0 0 1],"replicate","post"));
        My_pad = (padarray(My{i_t},[0 0 1],"replicate","post"));
        Mz_pad = (padarray(Mz{i_t},[0 0 1],"replicate","post"));
                
        %% compute the curl of magnetization (M)
        [curl_Mx,curl_My,curl_Mz,~] = curl(X,Y,Z,Mx_pad,My_pad,Mz_pad);
    
	    %% calculating vector potential (A) with concolution
        Ay_conv = A_*conv2(distance_kernel_c,curl_My(:,:,1),'valid');   % volume term
        Ay_conv = Ay_conv + A_*conv2(distance_kernel_b-distance_kernel_t,Mx{i_t},'valid')/dz;   % surface term
        Ay = Ay_conv(1:end-1,cpw_pos-cpw_pos(1)+1);

        AY(i_t,:,:) = Ay(:,:);
    end

    %% Approximate derivative with central differences
    dAY = diff(AY,1,1);

    dtime = diff(time);
    dtc = [dtime(1);dtime]+[dtime;dtime(end)];

    Ey = -([dAY(1,:,:,:);dAY]+[dAY;dAY(end,:,:,:)])./dtc;
   
    %% calculate the induced voltage in CPW
    %V1 = -sum(Ey(:,:,1),2)*1e-6;
    %V2 = -sum(Ey(:,:,2),2)*1e-6;
    
    V1 = -sum(Ey(:,:,1:10),[2,3])*(dy)/10;

    V2 = -sum(Ey(:,:,11:20),[2,3])*(dy)/10;

    V3 = -sum(Ey(:,:,21:30),[2,3])*(dy)/10;

    V = (V1+V3)/2-V2;

    %% fit sine curve
    sin_vec = sin(2*pi*f(i_f)*1e9*time);
    cos_vec = cos(2*pi*f(i_f)*1e9*time);
    re = sum(cos_vec.*V)*2/(length(time));
    im = sum(sin_vec.*V)*2/(length(time));

    V_real(i_f) = re;
    V_imag(i_f) = -im;

    %% plot signal and fit
    figure(i_f)
    plot(time*1e9,V1*1000,'-o');
    hold on
    plot(time*1e9,V2*1000,'-o');
    plot(time*1e9,V3*1000,'-o');
    plot(time*1e9,V*1000,'-o');
    plot(time*1e9,(re*cos_vec+im*sin_vec)*1000,'-k');
    hold off    
    grid on
    set(gca,'FontSize',16);
    xlabel('time [ns]')
    ylabel('deltaV [mV]')
    drawnow

end

%% plot resonance curves
figure(6000)
plot(f,V_real*1000,'.-','LineWidth',1)
hold on
plot(f,V_imag*1000,'.-','LineWidth',1)
%plot(f,abs(V_real+1i*V_imag)*1000,'.-','LineWidth',1)
hold off
grid on
xlabel('f (GHz)')
ylabel('Induced voltage (mV)')
set(gca,'FontSize',16,'LineWidth',1);
legend('real','imag')

%% plot resonance curves - with angle and abs
figure(6001)
plot(f,abs((V_real + 1i*V_imag)*1000),'.-','LineWidth',1)
grid on
xlabel('f (GHz)')
ylabel('abs(Z11) [mV]')
set(gca,'FontSize',16,'LineWidth',1);
legend('abs(Z11)')

figure(6002)
plot(f,angle((V_real + 1i*V_imag)*1000),'.-','LineWidth',1)
grid on
xlabel('f (GHz)')
ylabel('angle(Z11)')
set(gca,'FontSize',16,'LineWidth',1);
legend('angle(Z11)')

name = "voltage_vector_potential_Z11";
%saveas(gcf,dr+"/images/"+name+'.png')
%saveas(gcf,dr+"/images/"+name+'.fig')

save(dr+"/"+name+".mat","f","V_real","V_imag")
% % % save(dr+"/"+name+"_all"+".mat")
toc
profile report