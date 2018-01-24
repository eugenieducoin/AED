vdsat = 1.5;
vdlin = 0.1;
vdd = 1.5;
%start with transfer triode 
file1 = 'task1_IdVd_linear.csv';
format long
M = csvread(file1, 1, 0);

vg1 = M(:,1);
id_lin = M(:,2);

figure
hold on
plot(vg1, id_lin, 'k')

gm = diff(id_lin)./diff(vg1); %transconductance
plot(vg1(2:end), gm, 'b')

gm_max = max(gm)
gm_max_index = find(gm == gm_max);
vg_m = vg1(gm_max_index(1));

%max1 = find(vg1 == vg_m);
id_linm = id_lin(gm_max_index(1));

y_tangent = gm_max*vg1 - gm_max*vg_m + id_linm;
[min_tangent, index_min_tangent] = min(abs(y_tangent));
vth_lin = vg1(index_min_tangent)% + (vdlin/2)

plot(vg1, y_tangent, 'm')
title('Transfer characteristics - triode region')
legend('Id', 'Transconductance gm', 'Tangent to linear region of Id', 'Location', 'south')
hold off

%analysis for transfer saturation
file2 = 'task1_saturation.csv';
format long
D = csvread(file2, 1, 0);

vg2 = D(:,1);
id_sat = D(:,2);
figure
hold on
plot(vg2, id_sat, 'c')

id_satsq = sqrt(id_sat);
plot(vg2, id_satsq, 'm')

d = diff(id_satsq)./diff(vg2);
d_max = max(d);
d_max_ind = find(d==d_max);
vg_max2 = vg2(d_max_ind(1));

%i2ind = find(vg2 == vg_max2);
id_max2 = id_satsq(d_max_ind(1));

y = d_max*vg2 - d_max*vg_max2 + id_max2;
[miny2, indy2] = min(abs(y));
vthsat = vg2(indy2)

plot(vg2(1:end),y,'k')

title('Transfer characteristics - saturation region')
legend ('Id', 'Square Root of Id', 'Tangent to linear region of square root of id', 'Location', 'south')
hold off

%sub-threshold swing
log_id_sat = log10(id_sat);
N = length(vg2);
s_vec = zeros(1,(N-1));
delta_vg_vec = zeros(1,(N-1));
delta_log_id_vec = zeros(1,(N-1));

for t=1:(N-1)
    delta_vg_vec(t) = vg1(t+1)-vg1(t);
end    
for r=1:(N-1)
    delta_log_id_vec(r) = log10(id_lin(r+1)/id_lin(r));
end
for s1=1:(N-1)
    s_vec(s1) = delta_vg_vec(s1)/delta_log_id_vec(s1);
end

figure
hold on
semilogy(vg2, log_id_sat, 'b')
yyaxis right
plot(vg2(1:N-1), s_vec, 'g')
hold off
s = real(min(s_vec))*1000 %mV/dec

%on/off current ratio
a = vdd + vthlin - abs(vthlin);
b = vthsat - abs(vthsat);

[minA, indexMinA] = min(abs(vg2 - a));
%idon = id_lin(N);
idon = id_sat(indexMinA); 

[minB, indexMinB] = min(abs(vg2 - b));
idoff = id_sat(indexMinB);

on_off_current_ratio = idon/idoff

%dibl
num = abs(vth_lin-vthsat);
den = (vdsat-vdlin);
dibl_mV/V = num*1000/den

%analysis of output characteristics
file3 = 'task1_IdVd_sat.csv';
O = csvread(file3, 1, 0);

vd2 = O(:,1);
id_o = O(:,2);

figure
hold on
plot(vd2, id_o, 'k')

d = diff(id_o)./diff(vd2); %output conductance
[vd_sat, index_vd_sat] = min(abs(vd2 - vdsat));
go = d(index_vd_sat)
plot(vd2(2:end), d, 'b')
title('Output characteristics')
legend('Id', 'Location', 'south')
hold off
