%% Quick script for testing the TCP connection
fprintf(t,'Hello');
data = native2unicode(fread(t,13,'uchar'));
output = '';
for i = 1:length(data)-1
    output = [output data(i)];
end
data = string(output);
disp(data);