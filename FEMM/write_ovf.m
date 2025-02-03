function write_ovf(filename,title,dx,dy,dz,Nx,Ny,Nz,X,Y,Z,unit,label)

file = fopen(filename,'w');

fprintf(file,'# OOMMF OVF 2.0\n');
fprintf(file,'# Segment count: 1\n');
fprintf(file,'# Begin: Segment\n');
fprintf(file,'# Begin: Header\n');
fprintf(file,'# Title: %s GHz\n',title);
fprintf(file,'# meshtype: rectangular\n');
fprintf(file,'# meshunit: m\n');
fprintf(file,'# xbase: %E\n',dx/2);
fprintf(file,'# ybase: %E\n',dz/2);
fprintf(file,'# zbase: %E\n',dy/2);
fprintf(file,'# xstepsize: %E\n',dx);
fprintf(file,'# ystepsize: %E\n',dy);
fprintf(file,'# zstepsize: %E\n',dz);
fprintf(file,'# xnodes: %i\n',Nx);
fprintf(file,'# ynodes: %i\n',Ny);
fprintf(file,'# znodes: %i\n',Nz);
fprintf(file,'# xmin: 0\n');
fprintf(file,'# ymin: 0\n');
fprintf(file,'# zmin: 0\n');
fprintf(file,'# xmax: %E\n',Nx*dx);
fprintf(file,'# ymax: %E\n',Ny*dy);
fprintf(file,'# zmax: %E\n',Nz*dz);
fprintf(file,'# valuelabels: %s\n',[label,'_x ',label,'_y ',label,'_z']);
fprintf(file,'# valueunits: %s\n',[unit,' ',unit,' ',unit]);
fprintf(file,'# valuedim: %i\n',3);
fprintf(file,'# End: Header\n');
fprintf(file,'# Begin: Data Text\n');

fprintf(file,'%g %g %g\n',permute(cat(4,X,Y,Z),[4,2,1,3]));

fprintf(file,'# End: Data Text\n');
fprintf(file,'# End: Segment\n');

fclose(file);

