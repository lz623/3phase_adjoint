function [ ] = generate_well(schedule,nw)
gas_inject=schedule(:,3,:)==4;
schedule(:,2,:)=schedule(:,2,:)+177.1076*gas_inject.*schedule(:,2,:);
pro_t='../schedule.dat';
pro_t=fopen(pro_t,'w');
fprintf(pro_t,'%d\n',schedule(:,1,1));
fclose(pro_t);
schedules= '../wellschedule.dat';
constrain='../wellconstrain.dat';
control='../wellcontrol.dat';
schedules=fopen(schedules,'w');
constrain=fopen(constrain,'w');
control=fopen(control,'w');
for i=1:nw
    %schedule(:,i+1)=schedule(:,i+1)+0.1*schedule(:,i+1).*randn(nc);
    fprintf(schedules,'%d\n',schedule(:,2,i));
    
    fprintf(control,'%u\n',schedule(:,3,i));
    fprintf(constrain,'%u %u %u %u %u %u %u\n',schedule(:,4:10,i)');
    fprintf(schedules,'\n\n');
    fprintf(constrain,'\n\n');
    fprintf(control,'\n\n');
end
fclose(schedules);
fclose(control);
fclose(constrain);

end