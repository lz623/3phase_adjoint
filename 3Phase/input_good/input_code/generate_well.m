function [ ] = generate_well(well,schedule,nw)

wells= '../well_inform.dat';
wells=fopen(wells,'w');
fprintf(wells,'%d\n\n\n',nw);
for i=1:nw
    fprintf(wells,'%d\n',well(i,1));
    fprintf(wells,'%d\n',well(i,2));
    fprintf(wells,'%d\n',well(i,3));
    fprintf(wells,'%d\n',well(i,4));
    if well(i,1)
        fprintf(wells,'%u \t%u \t%u \t%u \t\n',well(i,5:8)');
    else
        fprintf(wells,'%f \t%f \t%f \t%f \t\n',well(i,5:8)');
    end
    fprintf(wells,'\n\n');
end
fclose(wells);
pro_t='../schedule.dat';
pro_t=fopen(pro_t,'w');
fprintf(pro_t,'%d\n',schedule(:,1));
fclose(pro_t);

schedules= '../wellschedule.dat';
constrain='../wellconstrain.dat';
control='../wellcontrol.dat';
schedules=fopen(schedules,'w');
constrain=fopen(constrain,'w');
control=fopen(control,'w');
for i=1:nw
    %schedule(:,i+1)=schedule(:,i+1)+0.1*schedule(:,i+1).*randn(nc);
    fprintf(schedules,'%d\n',schedule(:,(i-1)*9+2));
    fprintf(control,'%u\n',schedule(:,(i-1)*9+3));
    fprintf(constrain,'%u %u %u %u %u %u %u\n',schedule(:,(i-1)*9+4:(i-1)*9+10)');
    fprintf(schedules,'\n\n');
    fprintf(constrain,'\n\n');
    fprintf(control,'\n\n');
end
fclose(schedules);
fclose(control);
fclose(constrain);
end

