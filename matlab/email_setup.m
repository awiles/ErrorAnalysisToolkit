% email set up.

setpref('Internet','SMTP_Server','smtp.me.com');
setpref('Internet','E_mail','awiles@me.com');

%test
% clear msg;
% datetime = clock;
% filename = sprintf('freespace_%4d%02d%02d_%02d%02d', datetime(1:5));
% 
% msg = sprintf('Model computing complete.\n    Started:  %d-%02d-%02d %02d:%02d\n    Ended:    %d-%02d-%02d %02d:%02d\n    Filename: %s',...
%     datetime(1:5), datetime(1:5), filename);
% 
% sendmail('awiles@imaging.robarts.ca','FEMLAB Processing Complete',msg, 'readme.txt');