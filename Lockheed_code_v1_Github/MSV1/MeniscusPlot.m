h0 = 0.001;
R  = 0.250;
Lc = 0.020;

rmin = 0.0;
rmax = 1.0;
nr = 2001;
dr = (rmax-rmin)/(nr-1);

r = zeros(nr,1);
h = zeros(nr,1);
for i=1:nr
    r(i) = R*(rmin + (i-1)*dr);
    h(i) = h0*cosh(2*(R/Lc)*(r(i)/Lc))/cosh(2*R*R/Lc/Lc);
end

figure(10)
plot(r,h)