% Microstrip LNA Project Hand Calculations
% ECE4415
% Author: Curtis Mayberry
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

function Smith_project_example
clear all
close all
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('------')
s = [ (0.44*exp( j*137.1/180*pi)) (0.1 *exp( j* 9.5/180*pi)) ; ...
      (3.71*exp( j*16.7 /180*pi)) (0.29*exp(-j* 102.9/180*pi)) ];
fmin = 0.7; % dB
rOPT = (0.28*exp( j*179/180*pi));
rn   = 0.05;
  
% Check stability
[stable, k, del] = rollett_stability(s);
% del = calc_del(s);
% k =calc_k(s);
if(stable)
    disp(['Stability: Unconditionally Stable'])
else
    disp(['Stability: Not Unconditionally Stable']);
end
disp(['k = ' num2str(k) ', |del| = ' num2str(abs(del))])

% Setup planes
rS_plane_fig = figure(1);
% set(1,'Position',[1290 9 624 932]); % Communications Studio screen position
set(1,'Position',[-1911 9 944 988]); % Home screen position
rS_plane = gca;
smith;
title('\Gamma_S Plane: Minimum Noise Figure ','FontSize', 18)
rL_plane_fig = figure(2);
% set(2,'Position',[1929 9 624 932]); % Communications Studio screen position
set(2,'Position',[-951 9 944 988]); % Home screen position
rL_plane = gca;
smith;
title('\Gamma_L Plane: Minimum Noise Figure ','FontSize', 18)

% Stability Circles
% [out_centre, out_radius] = stability_circles(s, 'i',rS_plane)
% [in_centre, in_radius] = stability_circles(s, 'o',rL_plane)

% Maximum Stable Gain
MSG = 10*log10(abs(s(2,1))/abs(s(1,2)));
disp(['MSG = ' num2str(MSG) 'dB']);

%  Find Max Gp
Gp_max = (abs(s(2,1))/abs(s(1,2)))*(k-sqrt(k^2-1)); % Gonzalez eq 3.7.8/3.6.10
gp_max = Gp_max / (abs((s(2,1))))^2;
Gp_max_dB = 10*log10(Gp_max);
C2 = s(2,2)-del*conj(s(1,1));
rML = gp_max*conj(C2)/...
      (1+gp_max*((abs(s(2,2)))^2-(abs(del))^2));
  rMS = (calc_rin(s,rML))';
disp(['rML = ' num2str(abs(rML)) ' < ' num2str(angle(rML)*(180/pi))]);
disp(['Gp,max = Gt,max = Ga,max = ' num2str(Gp_max_dB) ' dB']);

% Gp Circles
%  Plot Constant Gp circles in rL plane
Gp_values = gp_circles(s,10*log10(maximum_gain(s)) - [1 2 3 ],rL_plane);
disp('Plotted Gp values (dB): ');
disp(Gp_values');

Gp_Label_Loc = [0.5+1i*0.33...
                0.58-1i*0.09...
                0.6-1i*0.38];
plot(rL_plane,rML,'m+')
text(real(rML)-.05,imag(rML)+.045, '\Gamma_M_L','Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(rML)-.05,imag(rML)-.045, ['G_P_,_m_a_x = ' num2str(Gp_max_dB,4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(1))),imag(z_to_r(Gp_Label_Loc(1))), ['G_P = ' num2str(Gp_values(1),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(2))),imag(z_to_r(Gp_Label_Loc(2))), ['G_P = ' num2str(Gp_values(2),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(3))),imag(z_to_r(Gp_Label_Loc(3))), ['G_P = ' num2str(Gp_values(3),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');

% Ga Circles
%  Plot Constant Ga circles in rS plane
Ga_values = ga_circles(s,10*log10(maximum_gain(s)) - [1 2 3 ],rS_plane);
disp('Plotted Ga values (dB): ');
disp(Ga_values');

Ga_Label_Loc = [0.48-1i*0.6...
                0.88-1i*0.66...
                1.4-1i*0.66];
plot(rS_plane,rMS,'m+')
text(real(rMS)-.05,imag(rMS)+.045, '\Gamma_M_S','Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(rMS)-.05,imag(rMS)-.045, ['G_A_,_m_a_x = ' num2str(Gp_max_dB,4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(1))),imag(z_to_r(Ga_Label_Loc(1))), ['G_P = ' num2str(Ga_values(1),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(2))),imag(z_to_r(Ga_Label_Loc(2))), ['G_P = ' num2str(Ga_values(2),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(3))),imag(z_to_r(Ga_Label_Loc(3))), ['G_P = ' num2str(Ga_values(3),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');

% Plot Constant NF circles in rS plane
noise_circles(fmin, rOPT, rn, (0.8:0.1:1) + 1e-5, rS_plane)
disp('Plotted NF values (dB): ');
disp((0.8:0.1:1)');

NF_Label_Loc = [0.46+1i*0.3...
                0.41+1i*0.4...
                0.373+1i*0.457];
text(real(rOPT)-.05,imag(rOPT)-.045, 'NF_m_i_n = 0.7dB','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(1))),imag(z_to_r(NF_Label_Loc(1))), 'NF = 0.8','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(2))),imag(z_to_r(NF_Label_Loc(2))), 'NF = 0.9','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(3))),imag(z_to_r(NF_Label_Loc(3))), 'NF = 1.0','Parent',rS_plane,'Color','blue','Interpreter','tex');




%%%%%%%%%%%%%%%%%
% See what happens with rS = rOPT, rL = rout' = rMLopt

rML_opt = conj(s(2,2)+(s(1,2)*s(2,1)*rOPT)/(1-s(1,1)*rOPT));
VSWRin_opt = calc_VSWRin(s,rOPT,rML_opt);
VSWRout_opt = calc_VSWRout(s,rOPT,rML_opt);
Gp_opt = 10*log10(calc_Gp(s,rML_opt));
Ga_opt = 10*log10(calc_Ga(s,rOPT));
Gt_opt = 10*log10(calc_Gp(s,rML_opt));
plot(rL_plane,rML_opt,'r+')
text(real(rML_opt)-.05,imag(rML_opt)+.045, '\Gamma_M_L_,_O_P_T','Parent',rL_plane,'Color','red','Interpreter','tex');


disp('--')
disp('Minimum Noise figure (rS = rOPT & rL = rML_opt)')
print_solution(s,rOPT,rML_opt,fmin,rn,rOPT)
% disp(['rOPT = ' num2str(abs(rOPT)) ' < ' num2str((180/pi)*angle(rOPT))]);
% disp(['rML_opt = ' num2str(abs(rML_opt)) ' < ' num2str((180/pi)*angle(rML_opt))]);
% disp(['VSWRin = ' num2str(VSWRin_opt)]);
% disp(['VSWRout = ' num2str(VSWRout_opt)]);
% disp(['Gp_opt = ' num2str(Gp_opt)]);

%%%%%%%%%%%%%%%%%
% See what happens with rS = rOPT, 
%  Now choose rL for a VSWR not equal to 1.
%  Then draw constant VSWRout = 2 circle and see if a rL on this circle
%   gives a reasonable VSWRin.

% VSWRout Circle
% Try VSWRout = 2
disp('--')
disp('Try VSWRout = 2')
% VSWRout Circle
VSWRout = 2;
[VSWRin_min, rL_VSWRin_min, theta_VSWRin_min] = plot_VSWRoutCircle(rL_plane,s,rOPT,VSWRout,8);
text(real(z_to_r(0.51+1i*0.7)),imag(z_to_r(0.51+1i*0.7)), '2.0','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
text(real(z_to_r(0.45+1i*0.55)),imag(z_to_r(0.45+1i*0.55)), 'VSWR_O_U_T','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
disp(['min VSWRin = ' num2str(VSWRin_min)]);
% Choose min. VSWRin
rL = rL_VSWRin_min;

plot(rS_plane, rOPT, 'r+');
text(real(rOPT)-.05,imag(rOPT)+.045, '\Gamma_O_P_T','Parent',rS_plane,'Color','red','Interpreter','tex');
plot(rL_plane, rL, 'r+');
text(real(rL)-.05,imag(rL)+.045, '\Gamma_L_,_O_P_T','Parent',rL_plane,'Color','red','Interpreter','tex');


disp('--')
disp('Choose theta = -pi/2')
disp(['rS = rOPT = ' num2str(abs(rOPT)) ' < ' num2str((180/pi)*angle(rOPT))]);
print_solution(s,rOPT,rL,fmin,rn,rOPT)

%%%%%%%%%%%%%%%%%
% Try a smaller value of VSWRout = 1.8
disp('--')
disp('Try VSWRout = 1.8')
% VSWRout Circle
VSWRout = 1.8;
[VSWRin_min, rL_VSWRout_min, theta_VSWRout_min] = plot_VSWRoutCircle(rL_plane,s,rOPT,VSWRout,8);
text(real(z_to_r(0.51+1i*0.765)),imag(z_to_r(0.51+1i*0.765)), '1.8','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
disp(['min VSWRin = ' num2str(VSWRin_min)]);

%%%%%%%%%%%%%%%%%
% Try a smaller value of VSWRout = 1.6
disp('--')
disp('Try VSWRout = 1.6')
% VSWRout Circle
VSWRout = 1.6;
[VSWRin_min, rL_VSWRout_min, theta_VSWRout_min] = plot_VSWRoutCircle(rL_plane,s,rOPT,VSWRout,8);
text(real(z_to_r(0.51+1i*0.83)),imag(z_to_r(0.51+1i*0.83)), '1.6','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
disp(['min VSWRin = ' num2str(VSWRin_min)]);

%%%%%%%%%%%%%%%%%%%%%%
% Final Solution
% Give up a little NF for a larger increase in gain
% Lose 0.1dB of NF, add 1dB of gain
% Trade off can be seen on rS plane as you move from rSopt
%  to rS
rS = rOPT - 0.1 - j*0.2;
% rS = rOPT - 0.088 - j*0.21; % Not really much of an improvement...

disp('  ')
disp('-----')
disp('Give up a little NF for a larger increase in gain')
disp('Lose 0.1dB of NF, add 1dB of gain')
disp(['Choose rS = ' num2str(abs(rS)) ' < ' num2str((180/pi)*angle(rS))])


% Setup a new set of plots
% Setup planes
rS_plane_fig = figure(3);
% set(3,'Position',[1290 9 624 932]); % Communications Studio screen position
set(3,'Position',[-1911 9 944 988]); % Home screen position
rS_plane = gca;
smith;
title('\Gamma_S Plane: Final Solution ','FontSize', 18)
rL_plane_fig = figure(4);
% set(4,'Position',[1929 9 624 932]); % Communications Studio screen position
set(4,'Position',[-951 9 944 988]); % Home screen position
rL_plane = gca;
smith;
title('\Gamma_L Plane: Final Solution ','FontSize', 18)

% Gp Circles
%  Plot Constant Gp circles in rL plane
Gp_values = gp_circles(s,10*log10(maximum_gain(s)) - [1 2 3 ],rL_plane);
disp('Plotted Gp values (dB): ');
disp(Gp_values);
Gp_Label_Loc = [0.54+1i*0.245...
                0.58-1i*0.09...
                0.6-1i*0.38];
plot(rL_plane,rML,'m+')
text(real(rML)-.05,imag(rML)+.045, '\Gamma_M_L','Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(rML)-.05,imag(rML)-.045, ['G_P_,_m_a_x = ' num2str(Gp_max_dB,4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(1))),imag(z_to_r(Gp_Label_Loc(1))), ['G_P = ' num2str(Gp_values(1),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(2))),imag(z_to_r(Gp_Label_Loc(2))), ['G_P = ' num2str(Gp_values(2),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Gp_Label_Loc(3))),imag(z_to_r(Gp_Label_Loc(3))), ['G_P = ' num2str(Gp_values(3),4)],'Parent',rL_plane,'Color','magenta','Interpreter','tex');

% Ga Circles
%  Plot Constant Ga circles in rS plane
Ga_values = ga_circles(s,10*log10(maximum_gain(s)) - [1 2 3 ],rS_plane);
disp('Plotted Ga values (dB): ');
disp(Ga_values);
Ga_Label_Loc = [0.48-1i*0.6...
                0.88-1i*0.66...
                1.4-1i*0.66];
plot(rS_plane,rMS,'m+')
text(real(rMS)-.05,imag(rMS)+.045, '\Gamma_M_S','Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(rMS)-.05,imag(rMS)-.045, ['G_A_,_m_a_x = ' num2str(Gp_max_dB,4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(1))),imag(z_to_r(Ga_Label_Loc(1))), ['G_P = ' num2str(Ga_values(1),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(2))),imag(z_to_r(Ga_Label_Loc(2))), ['G_P = ' num2str(Ga_values(2),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');
text(real(z_to_r(Ga_Label_Loc(3))),imag(z_to_r(Ga_Label_Loc(3))), ['G_P = ' num2str(Ga_values(3),4)],'Parent',rS_plane,'Color','magenta','Interpreter','tex');

% Plot Constant NF circles in rS plane
noise_circles(fmin, rOPT, rn, (0.8:0.1:1) + 1e-5, rS_plane)
disp('Plotted NF values (dB): ');
disp((0.8:0.1:1));
NF_Label_Loc = [0.46+1i*0.3...
                0.41+1i*0.4...
                0.373+1i*0.457];
plot(rS_plane,rOPT,'b+')
text(real(rOPT)-.05,imag(rOPT)-.045, 'NF_m_i_n = 0.7dB','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(1))),imag(z_to_r(NF_Label_Loc(1))), 'NF = 0.8','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(2))),imag(z_to_r(NF_Label_Loc(2))), 'NF = 0.9','Parent',rS_plane,'Color','blue','Interpreter','tex');
text(real(z_to_r(NF_Label_Loc(3))),imag(z_to_r(NF_Label_Loc(3))), 'NF = 1.0','Parent',rS_plane,'Color','blue','Interpreter','tex');


% VSWRout Circle
text(real(z_to_r(0.4+1i*0.53)),imag(z_to_r(0.4+1i*0.53)), 'VSWR_O_U_T','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
disp('--')
disp('Try VSWRout = 1.8')
VSWRout = 1.8;
[VSWRin_min, rL_VSWRout_min, theta_VSWRout_min] = plot_VSWRoutCircle(rL_plane,s,rS,VSWRout,8);
text(real(z_to_r(0.51+1i*0.86)),imag(z_to_r(0.51+1i*0.83)), '1.8','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
rML_rS = conj(s(2,2)+(s(1,2)*s(2,1)*rS)/(1-s(1,1)*rS));
plot(rL_plane,rML_rS, 'b+');
disp(['min VSWRin = ' num2str(VSWRin_min)]);
rL = rL_VSWRout_min;
plot(rS_plane,rS,'r+')

disp(' ')
disp('Solution:')
print_solution(s,rS,rL,fmin,rn,rOPT)

% VSWRout Circle
disp('--')
disp('Try VSWRout = 1.5')
VSWRout = 1.5;
[VSWRin_min, rL_VSWRout_min, theta_VSWRout_min] = plot_VSWRoutCircle(rL_plane,s,rS,VSWRout,16);
text(real(z_to_r(0.51+1i*0.76)),imag(z_to_r(0.51+1i*0.76)), '1.5','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');

disp(['min VSWRin = ' num2str(VSWRin_min)]);
rL = rL_VSWRout_min;
plot(rS_plane,rS,'r+')

rL_final_solution = rL;

disp(' ')
disp('Solution:')
print_solution(s,rS,rL,fmin,rn,rOPT)

% VSWRout Circle
disp('--')
disp('Try VSWRout = 1.3')
VSWRout = 1.3;
[VSWRin_min, rL_VSWRout_min, theta_VSWRout_min] = plot_VSWRoutCircle(rL_plane,s,rS,VSWRout,8);
text(real(z_to_r(0.51+1i*0.68)),imag(z_to_r(0.51+1i*0.68)), '1.3','Parent',rL_plane,'Color','blue','FontSize',7,'Interpreter','tex');
disp(['min VSWRin = ' num2str(VSWRin_min)]);
rL = rL_VSWRout_min;
plot(rS_plane,rS,'r+')

disp(' ')
disp('Solution:')
print_solution(s,rS,rL,fmin,rn,rOPT)

plot(rL_plane,rL_final_solution,'r*')

disp(' ')
disp('-----')
disp('FINAL SOLUTION:')
print_solution(s,rS,rL_final_solution,fmin,rn,rOPT)
text(real(rS)-.05,imag(rS)+.047, '\Gamma_S_,_f_i_n_a_l','Parent',rS_plane,'Color','red','Interpreter','tex');
text(real(rL_final_solution)-.05,imag(rL_final_solution)+.042, '\Gamma_L_,_f_i_n_a_l','Parent',rL_plane,'Color','red','Interpreter','tex');







%%%%%%%%%%%%%
% Functions %
%%%%%%%%%%%%%
function Gp_gammaL_to_gammaS(s,del,Cp_rl,Rp_rl,Gp)
Cp_rs = (((1-s(2,2)*Cp_rl)*(s(1,1)-del*Cp_rl)')-((Rp_rl)^2*(del')*s(2,2)))/...
        ((abs(1-(s(2,2)*Cp_rl)))^2-Rp_rl^2*abs(s(2,2))^2);
    
Rp_rs = (Rp_rl*abs(s(1,2)*s(2,1)))/...
        (abs(abs(1-s(2,2)*Cp_rl)^2-Rp_rl^2*abs(s(2,2))^2));

disp('--');
disp(['Gp = ' Gp 'dB']);
disp(['Cp_rs = ' num2str(abs(Cp_rs)) ' < ' num2str((180/pi)*angle(Cp_rs))]);
disp(['Rp_rs = ' num2str(abs(Rp_rs))]);

% Convert points in rL plane to pts in rS plane
function rS = rLplane_to_rSplane(rL,s)
rS = (s(1,1)+(s(1,2)*s(2,1)*rL)/(1-s(2,2)*rL))';

% Gonzalez eq. 3.8.6
function rb_mag = calc_rb_mag(rL,rout)
rb_mag = abs((rout-(rL'))/(1-rout*rL));

function del = calc_del(s)
del = s(1,1)*s(2,2) - S(1,2)*S(2,1);

function k = calc_k(s)
del = s(1,1)*s(2,2) - S(1,2)*S(2,1);
k = (1-(abs(s(1,1)))^2-(abs(s(2,2)))^2+(abs(del))^2)/...
    (2*abs(s(1,2)*s(2,1)));

function rin = calc_rin(s,rL)
rin = s(1,1)+(s(1,2)*s(2,1)*rL)/...
             (1-s(2,2)*rL);
         
% Gonzalez eq. 2.6.5
function rout = calc_rout(s,rS)
rout = s(2,2)+(s(1,2)*s(2,1)*rS)/...
              (1-s(1,1)*rS);  
function VSWRin = calc_VSWRin(s,rS,rL)
rin = calc_rin(s,rL);
ra_mag = abs((rin-rS')/(1-rin*rS));
VSWRin = (1+ra_mag)/(1-ra_mag);

function VSWRout = calc_VSWRout(s,rS,rL)
rout = calc_rout(s,rS);
rb_mag = abs((rout-rL')/(1-rout*rL));
VSWRout = (1+rb_mag)/(1-rb_mag);

function Gp = calc_Gp(s,rL)
rin = calc_rin(s,rL);
Gp = (1/(1-(abs(rin)).^2)*(abs(s(2,1))).^2*(1-(abs(rL)).^2)/((abs(1-s(2,2)*rL)).^2));

% Calculate ra/rb
function rab = calc_rab(VSWR)
rab = (VSWR-1)/(VSWR+1);

function Ga = calc_Ga(s,rS)
rout = calc_rout(s,rS);
Ga = ((1-abs(rS)^2)/(abs(1-s(1,1)*rS))^2)*(abs(s(2,1)))^2*(1/(1-(abs(rout))^2));

% Calculates the Noise Figure in dB
function F = calc_NF(Fmin, rn, rOPT, rS)
F = 10.^(Fmin/10) + (4*rn*(abs(rS-rOPT))^2)/...
           ((1-(abs(rS))^2)*(abs(1+rOPT))^2); % Gonzalez 4.3.4
 F = 10 * log10(F);
% yS = (1-rS)/(1+rS); % Gonzalez eq 4.3.2
% yOPT = (1-rOPT)/(1+rOPT); 
% F = 10.^(Fmin/10) + (rn/(real(yS)))*(abs(yS-yOPT))^2; % Gonzalez eq 4.3.4

function Gt = calc_Gt(s,rS,rL)
rin = calc_rin(s,rL);
Gt = ((1-abs(rS)^2)/abs(1-rin*rS)^2) * (abs(s(2,1))^2) * ((1-abs(rL)^2)/(abs(1-s(2,2)*rL)^2));

% function Gp_circles_on_rS_plane(s)

% Calculates the reflection coefficient of a normalized impedance z
function r = z_to_r(z)
r = (z-1)/(z+1);

% plots the VSWRoutCircle and selects the minimum VSWRin on that circle.
function [VSWRin_min, rL_VSWRin_min, theta_VSWRin_min] = plot_VSWRoutCircle(rL_plane,s,rS,VSWRout,N)
rb = calc_rab(VSWRout);
rout = calc_rout(s,rS);
Cvo = (rout'*(1-abs(rb)^2))/...
      (1-abs(rb*rout)^2);
Rvo = (abs(rb)*(1-abs(rout)^2))/...
      (1-abs(rb*rout)^2);
smith_circles(Cvo,Rvo,'b-',256,rL_plane)
% Test 8 points on the constant VSWRout circle
for idx = 1:N
    theta(idx) = -pi+(2*pi/N)*idx; % -pi+pi/4:pi/4:pi
    rL(idx) = Cvo+Rvo*exp(j*theta(idx));
    VSWRin(idx) = calc_VSWRin(s,rS,rL(idx));
    
    disp_theta{idx} = ['theta = ' num2str(theta(idx)/pi) ' pi rad'];
    disp_rL{idx} = ['rL = ' num2str(abs(rL(idx))) ' < ' num2str((180/pi)*angle(rL(idx)))];
    disp_VSWRin{idx} = ['VSWRin = ' num2str(VSWRin(idx))];
    disp_spacing{idx} = [32 32 32];
    plot(rL_plane, rL, 'b*');
%     disp('-')
%     disp(['theta = ' num2str(theta/pi) ' pi rad'])
%     disp(['rL = ' num2str(abs(rL)) ' < ' num2str((180/pi)*angle(rL))]);
%     disp(['VSWRin = ' num2str(VSWRin)])
end
% Display Table of pts on VSWRout circle
disp_theta = char(disp_theta);
disp_rL = char(disp_rL);
disp_VSWRin = char(disp_VSWRin);
disp_spacing = char(cellfun(@char,disp_spacing','UniformOutput',false));
disp([disp_VSWRin disp_spacing disp_rL disp_spacing disp_theta]);

% Calculate min VSWRin
[VSWRin_min min_idx] = min(VSWRin);
theta_VSWRin_min = theta(min_idx);
rL_VSWRin_min = rL(min_idx);
plot(rL_plane, rL_VSWRin_min, 'g*');

% Print Solution
% Prints important information about a given rS, rL solution
function print_solution(s,rS,rL,Fmin,rn,rOPT)
F = calc_NF(Fmin, rn, rOPT, rS);
VSWRin = calc_VSWRin(s,rS,rL);
VSWRout = calc_VSWRout(s,rS,rL);
Gp = calc_Gp(s,rL);
Ga = calc_Ga(s,rS);
Gt = calc_Gt(s,rS,rL);

disp(['rS = ' num2str(abs(rS)) ' < ' num2str((180/pi)*angle(rS))]);
disp(['rL = ' num2str(abs(rL)) ' < ' num2str((180/pi)*angle(rL))]);
disp(['VSWRin = ' num2str(VSWRin)])
disp(['VSWRout = ' num2str(VSWRout)])
disp(['NF = ' num2str(F) 'dB'])
disp(['Gp = ' num2str(10*log10(Gp)) 'dB'])
disp(['Ga = ' num2str(10*log10(Ga)) 'dB'])
disp(['Gt = ' num2str(10*log10(Gt)) 'dB'])

