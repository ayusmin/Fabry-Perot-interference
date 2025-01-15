n1 = 1; %real part of refractive index of medium 1
k1 = 0; %imaginary part of refractive index of medium 1
n2 = 2.4; %real part of refractive index of medium 2
k2 = 1.2; %imaginary part of refractive index of medium 2
n3 = 1; %real part of refractive index of medium 3
k3 = 0; %imaginary part of refractive index of medium 3
d = 2200; %separation between two films
a1 = 1; %absorption coefficient
H = 0.2.*d;
C = 0.08;
N1 = abs(n1 + k1*(0 + 1i)); %Complex refractive index of medium 1
N2 = abs(n2 + k2*(0 + 1i)); %Complex refractive index of medium 2
N3 = abs(n3 + k3*(0 + 1i)); %Complex refractive index of medium 3

theta1 = 0; %angle of incidence   
theta2 = abs(asin(N1*sin(theta1)/N2)); %angle of reflection from medium 1
theta3 = abs(asin(N2*sin(theta2)/N3)); %angle of reflection from medium 3


%reflection coefficients of s- polarized light for media 1,2 and 3
rs12 =((N1*cos(theta1))- (N2*cos(theta2)))/((N1*cos(theta1))+ (N2*cos(theta2)));
rs23 =((N2*cos(theta2))- (N3*cos(theta3)))/((N2*cos(theta2))+ (N3*cos(theta3)));
rs21 =((N2*cos(theta2))- (N1*cos(theta1)))/((N2*cos(theta2))+ (N1*cos(theta1)));

%reflection coefficients of p- polarized light for media 1,2 and 3
rp12 =((N1*cos(theta2))- (N2*cos(theta1)))/((N1*cos(theta2))+ (N2*cos(theta1)));
rp23 =((N2*cos(theta3))- (N3*cos(theta2)))/((N2*cos(theta3))+ (N3*cos(theta2)));
rp21 =((N2*cos(theta1))- (N1*cos(theta2)))/((N2*cos(theta1))+ (N1*cos(theta2)));

%transmission coefficients of s- polarized light for media 1,2 and 3
ts12 =(2*N1*cos(theta1))/((N1*cos(theta1))+ (N2*cos(theta2)));
ts23 =(2*N2*cos(theta2))/((N2*cos(theta2))+ (N3*cos(theta3)));
ts21 =(2*N2*cos(theta2))/((N2*cos(theta2))+ (N1*cos(theta1)));

%transmission coefficients of p- polarized light for media 1,2 and 3
tp12 =(2*N1*cos(theta2))/((N1*cos(theta2))+ (N2*cos(theta1)));
tp23 =(2*N2*cos(theta3))/((N2*cos(theta3))+ (N3*cos(theta2)));
tp21 =(2*N2*cos(theta1))/((N2*cos(theta1))+ (N1*cos(theta2)));

lambda = linspace(540, 900, 2103); %range of incident wavelength

%calculating the interefence function
for i = 1:length(lambda)
    IF = (n1^3.*(cos(theta1))^2.*C.*H.*...
        (1+(rs23).^2+(((2.*rs23)./((2.*pi.*n2.*H.*cos(theta2))./lambda(i))).*...
        cos((2.*pi.*n2.*cos(theta2).*(2.*d-H))./lambda(i)).*...
        sin((2.*pi.*n2.*H.*cos(theta2))./lambda(i)))))./...
        (2.*pi.*n2.*(n1.*cos(theta1)+n2.*cos(theta2)).^2.*...
        (1+(rs21.*rs23).^2-2.*rs21.*rs23.*cos((4.*pi.*n2.*d.*cos(theta2))./lambda(i))));
    IFs(i) = IF;

end
% IFs_norm = IFs./IFsn;
% If_abs = abs(IFs);

lambda_2 = (N2.*d)./lambda;
e_g = 1240./lambda; %bandgap of the materials

% IFs_norm_1 = IFs_norm - 0.7;
IFs_1 = IFs+0.5 ; %Interference fuction of s-polarized light

% IF_norm = (0.5.*(IFp_norm + IFs_norm))+0.7;

plot(lambda, IFs_1,'LineWidth',3);

xlim([540,900]);
ylim([0.1,1.05]);
xlabel('Wavelength (nm)');
ylabel('Interference Funtion (IF)');
title ('Fabry-Perot interference')



