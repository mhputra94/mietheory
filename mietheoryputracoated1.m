% Program Mie Theory untuk Partikel Bola Tunggal Coated
% Bagian Dalam Logam
% Bagian Luar Non Logam
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 22 Februari 2015

clear all
clc
% Input Material
matdal = 'silverpalik.dat'; % Material Dalam
nmatlu = 1.531; % Indeks Bias Material Luar
nling = 1.3334; % Indeks Bias Lingkungan
erei = 300:1:700; % Input Panjang Gelombang
rdal = 20; % jari - jari dalam
rlua = 50; % jari - jari luar
l = 50; % banyak iterasi
[ene, n, k] = textread(matdal, '%f %f %f','commentstyle','matlab');

units

enei = eV2nm ./ ene;
ni = spline( enei, n, erei );
ki = spline( enei, k, erei );
i = sqrt(-1);
nmat = ni + i*ki; % Indeks Bias Material Dalam

% Menghitung Mie Theory

x = (2*pi*nling*rdal)./erei;
y = (2*pi*nling*rlua)./erei;
m1 = nmat./nling;
m2 = nmatlu/nling;
k = (2*pi*nling)./erei;

o = length(erei);
sca = zeros(o,1);
ext = zeros(o,1);

% Menghitung Output
for j = 1:o
   for t = 1:l
       % Menentukkan Koefisien Hamburan
       
       % Menghitung Fungsi Riccati Bessel
      [psi, psit] = ricbestu(t,m2*x(j));
      [pso, psot] = ricbestu(t,m1(j)*x(j));
      [pse, pset] = ricbestu(t,y(j));
      [psa, psat] = ricbestu(t,m2*y(j));
      [cha, chat] = ricbesdua(t,m2*x(j));
      [che, chet] = ricbesdua(t,m2*y(j));
      [zet, zett] = ricbesga(t,y(j));
      chi = -1*cha; chit = -1*chat;
      cho = -1*che; chot = -1*chet;
      
      % Koefisien An => Ane dan Bn => Bne
      Anum = m2*psi*psot - m1(j)*psit*pso;
      Adenum = m2*chi*psot - m1(j)*chit*pso;
      Ane = Anum/Adenum;
      
      Bnum = m2*pso*psit - m1(j)*psi*psot;
      Bdenum = m2*chit*pso - m1(j)*psot*chi;
      Bne = Bnum/Bdenum;
      
      % Koefisien an dan bn
      anom = pse*(psat - (Ane*chot)) - ((m2*pset)*(psa - (Ane*cho)));
      adenom = zet*(psat - (Ane*chot)) - ((m2*zett)*(psa - (Ane*cho)));
      an = anom/adenom;
      
      bnom = ((m2*pse)*(psat - (Bne*chot))) - pset*(psa - (Bne*cho));
      bdenom = ((m2*zet)*(psat - (Bne*chot))) - zett*(psa - (Bne*cho));
      bn = bnom/bdenom;
      
      % Menghitung Hamburan dan Ekstinsi
      insca = (2*pi/(k(j))^2)*(2*t+1)*((abs(an))^2 + (abs(bn))^2);
       sca(j) = sca(j) + insca;
       inext = (2*pi/(k(j))^2)*(2*t+1)*((real(an + bn)));
       ext(j) = ext(j) + inext;
             
   end  
end

% Menghitung Serapan
abso = ext - sca;

% Plotting Grafik Hasil Mie Theory
plot(erei', ext,'r')
hold on
plot(erei', sca,'b')
plot(erei', abso,'g')
grid on
legend('Extinction','Scattering','Absorption')




