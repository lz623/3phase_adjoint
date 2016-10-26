%%Load input data
filename='reservoir_inform.xlsx';
sheet=1;
xlRange='B2:D2';
re_d = xlsread(filename,sheet,'B2:D2');
re_s=xlsread(filename,sheet,'B5:D5');
k_phi=xlsread(filename,sheet,'B19:B20');
sheet=2;
frac=xlsread(filename,sheet,'C2:D8');
sheet=3;
nw=xlsread(filename,sheet,'C1');
nsche=xlsread(filename,sheet,'B21');
well_range=['C4:J',num2str(3+nw)];
well=xlsread(filename,sheet,well_range);
idx=double('C');
sche_range=['C22:',char(idx+9*nw),num2str(nsche+21)];
schedule=xlsread(filename,sheet,sche_range);
%%Generate fracture system
%  generate_fracture(re_d(1,1),re_d(1,2),re_s(1,1),re_s(1,2),frac);
%%Generate perm poro field
generate_perm_por_field(re_d(1,1),re_d(1,2),k_phi(1),k_phi(2));
%%Generate well schedule
generate_well(well,schedule,nw)
%%Generate other input parameter
sheet=4;
n=xlsread(filename,sheet,'C1:C12');
parameter= ['../parameter.dat'];
parameter=fopen(parameter,'w');
fprintf(parameter,'%d\n',n);
fclose(parameter);