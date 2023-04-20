function total_reduction = doubleglass(f, para1, para2, d)

c_0 = 340;
rho_0 = 1.2;
f_c = 10^3 * (2.^((-3:2)));

r1 = para1.r;
r2 = para2.r;

m_total = para1.m + para2.m;
m_eff = para1.m*para2.m/(m_total);
K = rho_0*c_0^2/d;

f_0 = 1/(2*pi)*sqrt(K/m_eff); 

B_1 = para1.E *(para1.t^3/(12*(1-para1.nu^2)));
B_2 = para2.E *(para2.t^3/(12*(1-para2.nu^2)));
B_total = B_1 + B_2;

%f_c = c_0^2/(2*pi)*sqrt(m_total/B_total);
f_d = 55/d;

temp = zeros(1,length(f));
for i = 1:length(f)

    if f(i)<f_0 && f(i)<f_d
        temp(i) = 20*log10(2*pi*f(i).*m_total./(2*rho_0*c_0));

    elseif f(i)>=f_0 && f(i)<f_d
        temp(i) = r1(i) + r2(i) + 20*log10(f(i)*d) - 29;

    elseif f(i)>=f_0 && f(i)>=f_d
        temp(i) = r1(i) + r2(i) + 6;

    end

end

total_reduction = temp;

end