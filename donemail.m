function donemail(receiver_address,subject_string,message_string)
    my_default_email_address = 'glabmatemail@yahoo.com';
    SMTP_Username = 'glabmatemail';
    SMTP_Server = 'smtp.mail.yahoo.com';
    mypass = 'easypass123!@#';

    setpref('Internet','E_mail',my_default_email_address);
    setpref('Internet','SMTP_Server',SMTP_Server)
    setpref('Internet','SMTP_Username',SMTP_Username);
    setpref('Internet','SMTP_Password',mypass)

    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');

    sendmail(receiver_address,subject_string,message_string)
end