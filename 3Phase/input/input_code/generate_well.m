function [ ] = generate_well(well,nw)

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


end

