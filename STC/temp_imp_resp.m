function time_response = temp_Imp_Resp(n,k,t);

time_response = (k*t).^n.*exp(-k*t).*(1/factorial(n)-(k*t).^2/factorial(n+2));
