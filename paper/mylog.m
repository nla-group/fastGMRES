function mylog(s)

fid = fopen('mylog.txt', 'a+');
fprintf(fid, '%s', s);
fclose(fid);
fprintf(s);
