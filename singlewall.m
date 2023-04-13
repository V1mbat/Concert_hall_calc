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

temp = 20*log10(2*pi*f.*m./(2*rho_0*c_0)); 

for i = 1:length(f)
    if f(i) > f_c
        temp(i) = temp(i) + 10*log10(f(i)/f_c) + 10*log10(2*eta/pi);
    end
end
reduction_factor = temp;
end