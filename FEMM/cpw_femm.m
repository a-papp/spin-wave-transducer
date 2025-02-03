if ~exist('openfemm','file')
    addpath('C:\femm42\mfiles');
end

load("parameters.mat","dr")
[ ~, ~ ] = mkdir(dr);       % create folder if needed

f = 9.0;  % frequencies in GHz

%CPW geometries under investigation
frequencies =     f*1e9;
CPW_thicknesses = ones(size(frequencies))*0.16;
widths_ground =   ones(size(frequencies))*2;
widths_signal =   ones(size(frequencies))*2;
gaps =            ones(size(frequencies))*2.5;

for i_f = 1:length(frequencies)
    
    openfemm;
    newdocument(0);
    
    %% geometrical parameters in micrometers
    t_cpw = CPW_thicknesses(i_f);
    w_gnd = widths_ground(i_f);
    w_sgn = widths_signal(i_f);
    w_gap = gaps(i_f);
    
    %% enclosure dimensions in micrometers
    encH = 1000;
    encV = 1000;
    
    %% dimensions in mumax
    dx = 0.2e-6;
    dy = 1.0e-6;
    dz = 0.799e-6;
    
    Nx = 2500;
    Ny = 100;
    
    Xsize = dx*Nx;
    Ysize = dy*Ny;

    Ypos = -0.4;    % position of the centerline of YIG film
    
    %% define the striplines
    p_gnd = [0,        0
             w_gnd,    0
             w_gnd,    t_cpw
             0,        t_cpw];
    
    p_sgn = [w_gnd + w_gap,          0
             w_gnd + w_gap + w_sgn,  0
             w_gnd + w_gap + w_sgn,  t_cpw
             w_gnd + w_gap,          t_cpw];  

    p_gnd2 = [w_gnd + 2*w_gap + w_sgn,          0
              w_gnd + 2*w_gap + w_sgn + w_gnd,  0
              w_gnd + 2*w_gap + w_sgn + w_gnd,  t_cpw
              w_gnd + 2*w_gap + w_sgn,          t_cpw];  
    
    p_YIG = [-Xsize/2,  0
              Xsize/2,  0
              Xsize/2,  -dz
             -Xsize/2,  -dz]*1e6;
    
    DrawClosedPolygon(p_gnd);
    DrawClosedPolygon(p_sgn);
    DrawClosedPolygon(p_gnd2);
    DrawClosedPolygon(p_YIG);
    
    %% add materials
    mi_getmaterial('Copper')
    mi_modifymaterial('Copper', 5, 41); % Gold conductivity
    mi_getmaterial('Air');
        
    %% define circuits
    mi_addcircprop('Ground',5.0e-4,1);
    mi_addcircprop('Source',1.0e-3,1);
    
    %% add material labels
    mi_addblocklabel(w_gnd/2,t_cpw/2);
    mi_selectlabel(  w_gnd/2,t_cpw/2);
    mi_setblockprop('Copper', 0, dx*1e6, 'Ground', 0, 0, -1);
    mi_clearselected();
    
    mi_addblocklabel(w_gnd+w_gap+0.5*w_sgn,t_cpw/2);
    mi_selectlabel(  w_gnd+w_gap+0.5*w_sgn,t_cpw/2);
    mi_setblockprop('Copper', 0, dx*1e6, 'Source', 0, 0, 1);
    mi_clearselected();

    mi_addblocklabel(1.5*w_gnd + 2*w_gap + w_sgn,t_cpw/2);
    mi_selectlabel(  1.5*w_gnd + 2*w_gap + w_sgn,t_cpw/2);
    mi_setblockprop('Copper', 0, dx*1e6, 'Ground', 0, 0, -1);
    mi_clearselected();
    
    mi_addblocklabel(w_gap + 0.5*w_sgn + w_gnd,6.0*t_cpw);
    mi_selectlabel(w_gap + 0.5*w_sgn + w_gnd,6.0*t_cpw);
    mi_setblockprop('Air', 1, 0, 'None', 0, 0, 0);
    mi_clearselected();

    mi_addblocklabel(0,-t_cpw/2);
    mi_selectlabel(0,-t_cpw/2);
    mi_setblockprop('Air', 0, dx*1e6, 'None', 0, 0, 0);
    mi_clearselected();
        
    %% draw enclosure
    enc = [-encH/2, -encV/2
            encH/2, -encV/2
            encH/2,  encV/2
           -encH/2,  encV/2];
    DrawClosedPolygon(enc);
    
    enc = [enc;enc(1,:)];
    for i = 1:4
        enc_mid = (enc(i,:) + enc(i+1,:))/2;
        mi_selectsegment(enc_mid(1),enc_mid(2));
    end

    mi_addboundprop('Zero_A',0,0,0,0,0,0,0,0,0,0,0); % add boundary condition
    mi_setsegmentprop('Zero_A',0,0,0,0);
    mi_clearselected();
    
    %% solve
    mi_probdef(frequencies(i_f),'micrometers','planar',1e-8,100,30,0);
    mi_saveas(char(strrep(cd + "\CPW_design.fem", '\', '\\')));
    
    mi_analyze(0);
    mi_loadsolution();
    
    %% getting the flux density along YIG centerline
    B = mo_getb(((1:Nx)+0.5)*dx*1e6-250, Ypos*ones(1,Nx));

    %% generate the mumax input files for real and imag     
    filename_re = fullfile(dr,['CPW_DL_Real_',num2str(i_f),'.ohf']);
    write_ovf(filename_re, num2str(f(i_f)),dx,dy,dz, Nx,1,1, ...
              real(B(:,1)),zeros(Nx,1),real(B(:,2)), 'T','B' )   

    filename_im = fullfile(dr,['CPW_DL_Imag_',num2str(i_f),'.ohf']);
    write_ovf(filename_im, num2str(f(i_f)),dx,dy,dz, Nx,1,1, ...
              imag(B(:,1)),zeros(Nx,1),imag(B(:,2)), 'T','B' )


    closefemm;

end


%% function to draw polygons in FEMM
function DrawClosedPolygon(p)
    n = size(p,1);
    mi_addnode(p(1,1),p(1,2));
    for i = 2:n
        mi_addnode(p(i,1),p(i,2));
        mi_addsegment(p(i-1,1),p(i-1,2),p(i,1),p(i,2));
    end
    mi_addsegment(p(n,1),p(n,2),p(1,1),p(1,2));
end