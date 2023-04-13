function reduction_factor = singlewall(f, para)

t = para.t;
rho = para.rho;
eta = para.eta;
E = para.E;
nu = para.nu;
m = para.m;

c_0 = 340;
rho_0 = 1.2;

B = E *(t^3/(12*(1-nu^2)));

f_c = c_0^2/(2*pi)*sqrt(m/B);

R_0 = 20*log10(m) + 20*log10(f) - 42; 
perc = 0.0;
temp = zeros(1,length(f));

for i = 1:length(f)
    if f(i) < round((1-perc)*f_c)
        temp(i) = R_0(i) - 7;
    elseif f(i) >= round((1-perc)*f_c) && f(i) < round((1+perc)*f_c)
        temp(i) = R_0(i) - 7 + 10*log10(eta) + 8;
    elseif f(i) >= round((1+perc)*f_c)
        temp(i) = R_0(i) + 10*log10(f(i)/f_c - 1) + 10*log10(eta) - 2; 
    end
end
reduction_factor = temp;
end