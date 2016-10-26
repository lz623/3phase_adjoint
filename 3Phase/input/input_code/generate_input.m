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

for i=1:nw
    sche_range=['C',num2str((nsche+1)*(i-1)+22),':L',num2str((nsche+1)*i+20)];
     sche_range
    schedule(:,:,i)=xlsread(filename,sheet, sche_range);
end
%%Generate fracture system
%  generate_fracture(re_d(1,1),re_d(1,2),re_s(1,1),re_s(1,2),frac);
%%Generate perm poro field
generate_perm_por_field(re_d(1,1),re_d(1,2),k_phi(1),k_phi(2));
%%Generate well schedule
generate_well(well,nw)
generate_schedule(schedule,nw);
%%Generate other input parameter
sheet=4;
n=xlsread(filename,sheet,'C1:C12');
parameter= ['../parameter.dat'];
parameter=fopen(parameter,'w');
fprintf(parameter,'%d\n',n);
fclose(parameter);
% sheet=5;
% n_cs=xlsread(filename,sheet, 'B4')
% cs_range=['B8',':L',num2str(n_cs+8)];
% cs=xlsread(filename,sheet, sche_range);

