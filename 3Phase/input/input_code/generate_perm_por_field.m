%%This code is use to generate m_pri
function[m]=generate_perm_por_field(Nx,Ny,k,phi)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perm=k*ones(Nx*Ny,1);
perm=log(perm);
por=phi*ones(Nx*Ny,1);

m=[perm;perm;perm;por];
m_pri = ['../m_pri.dat'];
m_pri = fopen(m_pri,'w');
fprintf(m_pri,'%d\n',m');
fclose(m_pri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sche=linspace(stime,etime,30);
% sche=sche';
% timesche = ['schedule.dat'];
% fprintf(timesche,'%d\n',sche);
% fclose(timesche);

end