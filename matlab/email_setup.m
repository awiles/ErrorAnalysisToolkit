function email_setup(smtp, email)
% email set up.

setpref('Internet','SMTP_Server',smtp);
setpref('Internet','E_mail',email);