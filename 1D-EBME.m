% This is a parallel version of the 1-D Energy Balance Model modified to include the Greenhouse Effect for Earth-like planets. It finds and uses all of the available cores. (Compatible with the Slurm workload manager)

% This version calculates the fraction of the year the planet's surface is thawed at varying
% eccentricities, obliquities and SED as well as the edge of the Runway Green House (RGH)

% The code outputs results in an eccentricity vs. installation array in the .csv format 

clear all
feature('numcores')

% Flag for the the XUV flux calculate and mass loss (1 = ON, 0 = OFF)
%calcXUV = 1;

%Define Orbital Dynamics Parameters
obl = 0 %[0 45 90] 		% Define Obliquity Array
eccen = [0] %[0 0.3 0.8] 		% Define Eccentricity Array
ins=1:1:201; 			% Define the instellation array as percent of the Solar Constant i.e. ins = 100 = 100% S_0
per = 102.0651294475608;		% Argument of Periastron in degrees

% Define RunLength and Temporal Resolution
runlength=50;	% in years
ns=360;		% Number of time steps per year in days

% Radiative Transfer Parameters
A=203.3;
B=2.09;
Dmag = 0.44;
Toffset = -30;
p=num2str(per);

hadleyflag = 1.;		% Turn on Hadley Circulation
albedoflag = 0.;		% Make climatological albedo if no albedo feedback
ice_model = 1.;		% Turn on the Explicit Ice Model
coldstartflag = 0.; 	% Initial the Planet with Temperatures offset by Toffset

% Heat Capacities and Land-Ocean Coupling Constant
Cl = 0.45;
Cw = 9.8;
nu = 3;

% Output variables
P=zeros(length(eccen),length(ins));
TOLT = P;
MDOT = P;
RGHplanets=P;
Hind=P;

% Time loop parameters.
ts=90; % first day of run is present day equinox
tf=runlength;
nstepinyear=ns;
delt=1./nstepinyear;
nts=floor(tf/delt);
n_out=1:nts;

%size of domain.
jmx=50;	% Latitudinal resolution
jmx=2*floor(jmx/2);

%Parameters for mass loss fie to XUV flux
Rp = 6.371e6; 	%Planet's Radius
G = 6.674e-11;
Mp = 5.987e24;	%Planet's Mass
Ktide = 1;
ep = 0.1;
TO = 1.39e21; 	%Terrestrial ocean mass 

% Land Fraction
fl=0.01*ones(jmx,1);

% Ocean fraction.
fw = 1-fl;

nu_fw = nu./fw;
nu_fl = nu./fl;

%Assign Albedos for chosen stellar type
for rr = [2]
if rr == 1	% M-dwarf
    Aw=0.234;
    Al=0.332;
    Ai=0.315;
    zz = 'M';
end
if rr == 3	% K-dwarf
    Aw=0.302;
    Al=0.401;
    Ai=0.477;
    zz = 'K';
end
if rr == 2	% G-dwarf
    Aw=0.319;
    Al=0.415;
    Ai=0.514;
    zz = 'G';
end
if rr == 4	% F-dwarf
    Aw=0.329;
    Al=0.414;
    Ai=0.537;
    zz = 'F';
end

% freezing temperature of ocean in deg C
% 83 is ratio of Lf/Cp and 50 is depth of ML
Tfrz=-2;
conduct=2;
Lfice=9.8*83.5/50; % this is about right

%set up x array.
delx = 2.0/jmx;
x = [-1.0+delx:delx:1.0-delx]';
xfull=(-1+delx/2:delx:1-delx/2)';
phi = asin(xfull)*180/pi;

%heat diffusion coefficient.
if (hadleyflag)
  D=Dmag*(1+9*exp(-(x/sin(25*pi/180)).^6));
else
  D=Dmag*ones(length(x),1);
end

Cw_delt=Cw/delt;
Cl_delt=Cl/delt;
delt_Lf=delt/Lfice;

lambda=D/delx/delx.*(1-x.*x);
a=[0; -lambda];
c=[-lambda; 0];
b=-a-c;
Diff_Op =-( diag(b) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1));

tic

[d1,d2] = size(ins);
    z = num2str(obl);
    if obl == 0
       calcXUV = 1;  %XUV calculation is default for obl = 0 !
    else
       calcXUV = 0;  %XUV and mass loss calculation code incompatible with nonzero obliquities
    end

for ind2=ins
    scaleQ=(ind2-1)/200;

% Open files to write output data
fid1 = fopen(['Mdot_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid2 = fopen(['TOlt_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid3 = fopen(['RGHplanets_' zz '_' z '_aq_' p 'aq.csv'],'wb');
fid4 = fopen(['thaw_frac_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid5 = fopen(['habi_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid6 = fopen(['meantempglob_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');

for o=eccen
    ecc=(o-1)/length(eccen);

% Setup Initial Conditions
r=zeros(2*jmx,1);

bw=Cw_delt+B+nu_fw-a-c;
bl=Cl_delt+B+nu_fl-a-c;

Mw = diag(bw) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1);
Ml = diag(bl) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1);

M=zeros(2*jmx,2*jmx);
for j=1:jmx
  M(2*j-1,1:2:2*jmx)=Ml(j,:);
  M(2*j-1,2*j)=-nu_fl(j);
  M(2*j,2:2:2*jmx)  =Mw(j,:);
  M(2*j,2*j-1)=-nu_fw(j);
end

invM=inv(M);

L_out = zeros(jmx,length(n_out));
W_out = L_out;
if(ice_model), h_out=W_out; end

%set up inital T profile
if (coldstartflag), Toffset=-40; else, Toffset = 0.; end
L = 7.5+20*(1-2*xfull.^2)+Toffset;
W = 7.5+20*(1-2*xfull.^2)+Toffset;

%Initial Ice Cover
ice=find(W<Tfrz);
notice=find(W>=Tfrz);
h=zeros(jmx,1); k=h;
h(ice)=2;

% make climatological albedo if no albedo feedback
if (albedoflag)
  load temperatures.mat 
  clim_alb_l=zeros(jmx,360);
  clim_alb_w=clim_alb_l;
  n=1;
  for t = thedays
    [clim_alb_l(:,t),clim_alb_w(:,t)]=albedo_seasonal(Lann(:,n),Wann(:,n),xfull);
    n=n+1;
  end
end

%Set up initial albedos
if (albedoflag)
   alb_l=clim_alb_l(:,thedays(1)); alb_w=clim_alb_w(:,thedays(1));
else
   [alb_l,alb_w]=albedo_seasonal(L,W,xfull,Aw,Al,Ai);
end

%obtain annual array of daily averaged-insolation.
[G_xuv_F,M_xuv_F,insol] = seasonal_solar_XUV(xfull,obl,ecc,per);
insol=scaleQ*insol(:,[360 1:359]);

if (ice_model)
% compute a heat flux from the ocean 
% the flux asymtotes to 2 W/m2 as the ice area
% decreases down to a single gridbox of area (delx)
% and asymtotes to 0 as the ice area covers the hemisphere
% the ocean temperature is adjusted to
  S=insol(:,ts);
  rprimew=A-(1-alb_w).*S;
  rprimel=A-(1-alb_l).*S;
  r(1:2:2*jmx,1)=L*Cl_delt-rprimel;
  r(2:2:2*jmx,1)=W*Cw_delt-rprimew;
     k=zeros(jmx,1); Fnet=k;
   k(ice)=conduct./h(ice);
   r(2*ice)=k(ice)*Tfrz-rprimew(ice);
   dev=zeros(2*jmx,1);
   dev(2*ice)=-Cw_delt+k(ice);
   Mt=M+diag(dev);
   I=inv(Mt)*r;
   T=I;
   T(2*ice)=min(Tfrz,I(2*ice));
   L=T(1:2:2*jmx);
   W=T(2:2:2*jmx);
   I=I(2:2:2*jmx);
   Fnet(ice)=Diff_Op(ice,:)*I-rprimew(ice)-B*W(ice)...
       -nu_fw(ice).*(I(ice)-L(ice));Fnet=Fnet(:);

nhice=ice(find(ice>jmx/2+1));         % nh grid cells with ice
shice=ice(find(ice<jmx/2));
nhocn=notice(find(notice>jmx/2+1));   % nh grid cells with ocn
shocn=notice(find(notice<jmx/2));
nhicearea=sum(fw(nhice));             % nh ice area
shicearea=sum(fw(shice));
nhmax=sum(fw(jmx/2+1:jmx));               % nh max possible ocn/ice area
shmax=sum(fw(1:jmx/2));
nhfw=2*min(2-2*(nhicearea-delx)/nhmax,2);   % nh fw under ice
shfw=2*min(2-2*(shicearea-delx)/shmax,2);
Fnet(nhice)=Fnet(nhice)+nhfw;
Fnet(shice)=Fnet(shice)+shfw;

nhocnarea=nhmax-nhicearea;                % nh ice-free ocn area
shocnarea=shmax-shicearea;
nhdW=nhfw*nhicearea./nhocnarea/Cw_delt;   % nh adjust water temp
shdW=shfw*shicearea./shocnarea/Cw_delt;
W(nhocn)=W(nhocn)-nhdW;
W(shocn)=W(shocn)-shdW;

else
  ice=[];
end

idx_out=1;
tday=zeros(ts +2 +1+(nts-1)*360*delt,1);
yr=zeros(floor((-1+ts +2 +1+(nts-1)*360*delt)/360),1);
day=zeros(floor(ts +2 +1+(nts-1)*360*delt- floor((-1+ts +2 +1+(nts-1)*360*delt)/360)*360),1);
alb_out = zeros(jmx,nts);

%%%%% Begin main loop %%%%

for n = 1:nts
     tday(n)=ts +2 +1+(n-1)*360*delt;
     yr(n)=floor((-1+tday(n))/360);
     day(n)=floor(tday(n)- yr(n)*360);
     
%create initial albedo.

  if (albedoflag)
    nn = day(n);
    nn=nn-360*delt; if (nn<0), nn=360-360*delt *0.5; end
    alb_l=clim_alb_l(:,nn); alb_w=clim_alb_w(:,nn);
  else
    [alb_l,alb_w]=albedo_seasonal(L,W,xfull,Aw,Al,Ai);
  end

%Calculate insolation
   S = insol(:,day(n));
   gh = find(W > 46.2);   % find moist locations
   %rgh = find(W > 374); % find runaway locations
%Source terms
   rprimew=A-(1-alb_w).*S;
   rprimel=A-(1-alb_l).*S;
   %%insert!!%%
   A1 = 300;
   %B1 = 5;
   %A2 = 300-(B1*374);
   rprimew(gh)=A1-(1-alb_w(gh)).*S(gh);
   rprimel(gh)=A1-(1-alb_l(gh)).*S(gh);
   %rprimew(rgh)=A2-(1-alb_w(rgh)).*S(rgh);
   %rprimel(rgh)=A2-(1-alb_l(rgh)).*S(rgh);
   %%insert!!%%
   r(1:2:2*jmx,1)=L*Cl_delt-rprimel;
   r(2:2:2*jmx,1)=W*Cw_delt-rprimew;
  
  if (ice_model)
  % first consider where sea ice already exists
  ice=find(h>0.001);
  notice=find(h<=0.001);
  %if (length(ice)>0),
   k=zeros(jmx,1); Fnet=k;
   k(ice)=conduct./h(ice);
   r(2*ice)=k(ice)*Tfrz-rprimew(ice);
   dev=zeros(2*jmx,1);
   dev(2*ice)=-Cw_delt+k(ice);
   %%instert!!%%
   bw=Cw_delt+B+nu_fw-a-c;
   bl=Cl_delt+B+nu_fl-a-c;
   bw(gh)=Cw_delt+nu_fw(gh)-a(gh)-c(gh);
   bl(gh)=Cl_delt+nu_fl(gh)-a(gh)-c(gh);
   %bw(rgh)=Cw_delt+B1+nu_fw(rgh)-a(rgh)-c(rgh);
   %bl(rgh)=Cl_delt+B1+nu_fl(rgh)-a(rgh)-c(rgh);

   Mw = diag(bw) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1);
   Ml = diag(bl) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1);

   M=zeros(2*jmx,2*jmx);
   for j=1:jmx
     M(2*j-1,1:2:2*jmx)=Ml(j,:);
     M(2*j-1,2*j)=-nu_fl(j);
     M(2*j,2:2:2*jmx)  =Mw(j,:);
     M(2*j,2*j-1)=-nu_fw(j);
   end
    %%insert!!%%
   Mt=M+diag(dev);
   I=inv(Mt)*r;
   T=I;
   T(2*ice)=min(Tfrz,I(2*ice));
   L=T(1:2:2*jmx);
   W=T(2:2:2*jmx);
   I=I(2:2:2*jmx);
   Fnet(ice)=Diff_Op(ice,:)*I-rprimew(ice)-B*W(ice)...
       -nu_fw(ice).*(I(ice)-L(ice));Fnet=Fnet(:);

% compute a heat flux from the ocean

% the flux asymtotes to 2 W/m2 as the ice area
% decreases down to a single gridbox of area (delx)
% and asymtotes to 0 as the ice area covers the hemisphere
% the ocean temperature is adjusted to

nhice=ice(find(ice>jmx/2+1));         % nh grid cells with ice
shice=ice(find(ice<jmx/2));
nhocn=notice(find(notice>jmx/2+1));   % nh grid cells with ocn
shocn=notice(find(notice<jmx/2));
nhicearea=sum(fw(nhice));             % nh ice area
shicearea=sum(fw(shice));
nhmax=sum(fw(jmx/2+1:jmx));               % nh max possible ocn/ice area
shmax=sum(fw(1:jmx/2));
nhfw=2*min(2-2*(nhicearea-delx)/nhmax,2);   % nh fw under ice
shfw=2*min(2-2*(shicearea-delx)/shmax,2);
Fnet(nhice)=Fnet(nhice)+nhfw;
Fnet(shice)=Fnet(shice)+shfw;

nhocnarea=nhmax-nhicearea;                % nh ice-free ocn area
shocnarea=shmax-shicearea;
nhdW=nhfw*nhicearea./nhocnarea/Cw_delt;   % nh adjust water temp
shdW=shfw*shicearea./shocnarea/Cw_delt;
W(nhocn)=W(nhocn)-nhdW;
W(shocn)=W(shocn)-shdW;
   h(ice)=h(ice)-delt_Lf*Fnet(ice);
   h(ice)=max(0.,h(ice));
%  end

% second consider the conditions of new ice growth over the ocean
  cold=find( T(2*notice)<Tfrz);
  new=notice(cold);
  if (length(new)>0)
    h(new)=-Cw/Lfice*( W(new)-Tfrz );
    W(new)=Tfrz;
  end

 else
   T=invM*r;
   L=T(1:2:2*jmx);
   W=T(2:2:2*jmx);
 end

%output
alb_out(:,n)= fw.*alb_w(:)+fl.*alb_l(:);
   if n==n_out(idx_out) 
     L_out(:,idx_out) = L(:);
     W_out(:,idx_out) = W(:);
     if(ice_model), h_out(:,idx_out) = h(:); end
     idx_out=idx_out+1;
   end
end

n=idx_out-1;
days=(nstepinyear-1):-1:0;
tt=n-days;
thedays=day(tt);
Lann=L_out(:,tt);
Wann=W_out(:,tt);

MGH = (Wann > 47 & Wann < 77);   %Moist Greenhouse
MGHL = (Wann > 77 & Wann < 374); %Running with large H2O mixing ratios
%RGH = (Wann > 374);              %Late Runaway Greenhouse
h_lastyear = h_out(:,(10441:10800));
alb_lastyear = alb_out(:,(10441:10800));

% Determine if climate is in a runaway.
Tg50=(L_out(:,(17641:18000))'*fl+W_out(:,(17641:18000))'*fw)/jmx;
Tg40=(L_out(:,(14041:14400))'*fl+W_out(:,(14041:14400))'*fw)/jmx;

if (Tg50 - Tg40) < 1 % Planet in a RGH
    eq = 1;
else
    eq = 0;
    RGHplanets(o,ind2) = 1;
end
if Tg50 > 370
    RGHplanets(o,ind2) = 1;
end

%Calculate thaw fraction / habitability index
m=0;
H=zeros(50,360);
Hab=zeros(50,360);
Mdot = zeros(ns);

for i=1:1:360
    K=0;
    for j=1:1:50
        if h_lastyear(j,i)<=10^(-2)
            K=K+1;
        end
        if 0 < Wann(j,i) < 100
            H(j,i)=1;
        end
    end
    Hab(:,i) = H(:,i) .* cos(phi*pi/180);
    if K>=17
        m=m+1;
    end
    if calcXUV == 1
    %Calculate the mass loss while in the MGH/RGH state due to the XUV flux
    %RGHlat = RGH(:,i);
    MGHLlat = MGHL(:,i);
    MGHlat = MGH(:,i);
    if zz == 'G'
    XUVday = G_xuv_F(:,i) * scaleQ;
    elseif zz == 'M'
    XUVday = M_xuv_F(:,i) * scaleQ;
    end
    %Find the are and flux in the RGH
    %lat3 = find(RGHlat);
    lat2 = find(MGHlat);
    lat = find(MGHLlat);
    l2 = isempty(lat2);
    l3 = isempty(lat);
    %l1 = isemtpy(lat3);
    [dx2,dy2]=size(lat2);      % Find the normalization constant ~ the number of latitude entries in MGH state
    if dx2~=1     
        jj2=dx2;
    else
        jj2=dy2;
    end
    [dx,dy]=size(lat);   % Find the normalization constant ~ the number of latitude entries in MGHL state
    if dx~=1     
        jj=dx;
    else
        jj=dy;
    end
    XUVV1 = XUVday(:);
    XUVV2 = XUVday(MGHlat);
    XUVV3 = XUVday(MGHLlat);
    %XUVV4 = XUVday(RGHlat);
    XUVsumMGH = sum(XUVV2(:)); % Total XUV flux absorved by MGH lats
    XUVsumMGHL = sum(XUVV3(:));
    XUVsumRGH = sum(XUVV1(:));
    %XUVsumRGH2 = sum(XUVV4(:));

    if eq == 0 %Planet during the runaway
        Mdot(i) = ep .* pi * Rp^3.* XUVsumRGH./jmx ./ (G .* Mp .* Ktide);
    elseif (Tg50 > 370) %Planet in RGH
        Mdot(i) = ep .* pi * Rp^3.* XUVsumRGH./jmx ./ (G .* Mp .* Ktide);
    elseif ((l3 == 1) & (l2 == 1)) % If neither is present
        ph  = 0;
        ph2 = 0;
        Mdot(i) = 0;
    elseif ((l3 == 1) & (l2 == 0)) % if only MGH
        ph  = 0;
        la2 = min(lat2);
        maxphiMGH = abs(phi(la2));
        ph2 = maxphiMGH.*pi./180;
        Mdot(i) = ep/1000 .* pi *Rp^3*((2*ph2 + 2*sin(pi/2-ph2) *cos(pi/2-ph2)) - (2*ph + 2*sin(pi/2-ph) *cos(pi/2-ph)))...
        .* XUVsumMGH./jj2 ./ (G .* Mp .* Ktide);
    elseif ((l3 == 0) & (l2 == 0)) % IF both present
        la = min(lat);
        maxphiMGHL = abs(phi(la));  % Top edge of the MGHL zone
        ph = maxphiMGHL.*pi./180;
        la2 = min(lat2);
        maxphiMGH = abs(phi(la2));
        ph2 = maxphiMGH.*pi./180;  % Top edge of the MGH zone
        assert(ph2 > ph);
        Mdot(i) = ep .* pi * Rp^3*(2*ph + 2*sin(pi/2-ph) *cos(pi/2-ph)).* XUVsumMGHL./jj ./ (G .* Mp .* Ktide)...
        + ep/1000 .* pi *Rp^3*((2*ph2 + 2*sin(pi/2-ph2) *cos(pi/2-ph2)) - (2*ph + 2*sin(pi/2-ph) *cos(pi/2-ph)))...
        .* XUVsumMGH./jj2 ./ (G .* Mp .* Ktide);
    elseif ((l3 == 0) & (l2 == 1)) %Only MGHL
        ph2 = 0;
        la = min(lat);
        maxphiMGHL = abs(phi(la));  % Top edge of the RGH zone
        ph = maxphiMGHL.*pi./180;
        Mdot(i) = ep .* pi * Rp^3*(2*ph + 2*sin(pi/2-ph) *cos(pi/2-ph)).* XUVsumMGHL./jj ./ (G .* Mp .* Ktide);
    end
    end
end

Mlosstot = sum(Mdot(:)) .* 3.154e7; % kg lost per year
TOlt = TO ./ Mlosstot ./ 1.0e6; % time to lose entire ocean in Myr
MDOT(o) = Mlosstot;
TOLT(o) = TOlt;

P(o,ind2) = m/360;
sd = cos(phi*pi/180);
norm = sum(sd(:))*360;
Hind(o) = sum(Hab,'all')/(norm);

% Shortwave absorbed
%asr = (1-alb_lastyear).*insol; % (360,jmx=50)
%ASR = asr > 282;
%GH = find(ASR);
%[GHw,wh] = size(GH);
meantempglob(o) = mean2(Wann);
%disp([num2str(o) ',' num2str(ind2-1)]);
end

fprintf(fid1,'%f\t',Mdot,ind2);
fprintf(fid1,'\n');
fprintf(fid2,'%f\t',TOlt,ind2);
fprintf(fid2,'\n');
fprintf(fid3,'%f\t',RGHplanets,ind2);
fprintf(fid3,'\n');
fprintf(fid4,'%f\t',P,ind2);
fprintf(fid4,'\n');
fprintf(fid5,'%f\t',Hind,ind2);
fprintf(fid5,'\n');
fprintf(fid6,'%f\t',meantempglob,ind2);
fprintf(fid6,'\n');
fid1 = fopen(['Mdot_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid2 = fopen(['TOlt_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid3 = fopen(['RGHplanets_' zz '_' z '_aq_' p 'aq.csv'],'wb');
fid4 = fopen(['thaw_frac_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid5 = fopen(['habi_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
fid6 = fopen(['meantempglob_RGH_' zz '_' z '_ao_' p 'aq.csv'],'wb');
disp([num2str((ind2-1))]);
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
end
toc
%if calcXUV == 1
%mdot = cell2mat(MDOT);
%tolt = cell2mat(TOLT);
%RGHplan = cell2mat(RGHplanets);
%csvwrite(['Mdot_RGH_' zz '_' z '_ao_' p 'aq.csv'],mdot);
%csvwrite(['TOlt_RGH_' zz '_' z '_ao_' p 'aq.csv'],tolt);
%csvwrite(['RGHplanets_' zz '_' z '_aq_' p 'aq.csv'],RGHplan);