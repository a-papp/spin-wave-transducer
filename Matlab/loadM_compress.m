N=1000; % number of timesteps to load
load("parameters.mat","f","dr")
time=zeros(N,1);
[Mx{N}, My{N}, Mz{N}] = deal(0);
for i_f = 1:length(f)
    for i = 1:N
        num=num2str(i);
        while(size(num,2)<6);num=strcat('0',num);end
        mfile=strcat(dr,'/mumax_f',num2str(i_f),'.out/','m_xrange1000-1500_zrange0_',num,'.ovf');

        data = read_ovf(mfile);
        time(i) = data.time;
    
        Mx{i}=data.X;
        My{i}=data.Y;
        Mz{i}=data.Z;
    end    
    dx = data.dx;
    dy = data.dy;
    dz = data.dz;
    I = data.I;
    J = data.J;
    K = data.K;
    save([dr,'/results_compressed_f',num2str(i_f)],...
          "Mx","My","Mz","dx","dy","dz","I","J","K","time")
end
