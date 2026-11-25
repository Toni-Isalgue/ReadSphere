%
% This is a development program. From a matrix of float values
% (Temperatures)
% coming from a spherical thermography, produces:
% plot (false color picture), 
% violin plots (lineal and cos-weighted), (density distribution of points),
% averages and standard deviations, and
% decomposition in spherical harmonics (useful to obtain the 
% preferred orientation of radiation)

% Input: Data as in a csv file renamed to Dataexample.dat, 
% using decimal point for the numbers and comma as separators

% example data - Contains temperature as function of theta (different
% files= different theta value, starts at zenit), and fi (different rows,
% from starting fi=0 at left of the "picture"), theta: 0 -> pi rad; fi: 0
% -> 2 pi rad

A1=readmatrix('f:\ReadSphere\Dataexample.dat');
pid2=pi/2;

% for the violin plots, re-dimension
%  and avoid the nan

A= reshape(A1, [], 1);
i1=isnan(A);
A(i1)=[];

% figure foto from Dataexample.dat:
% Eventually, ask for temperature limits for the tretment 
% minI1=input("minI1?");
% maxI1=input("maxI1?");
% Values indicated for Dataexample.dat:

minI1=-35;
maxI1=85;
figure(1); %figure of first thermo pictures, "as is"
clims=[minI1 maxI1];
imagesc(A1, clims);
colorbar; %Scale of temperatures

% Take into account solid angle covered by dat: cos weight. Attention!

% Open/Create results file

fiop=fopen('f:\ReadSphere\Results.dat','w');
formatesp='No. elements, lineal average, lineal stdev: %d %d %5.3f %5.3f \n';
fprintf(fiop,'Data in Dataexample.dat \n');
fprintf(fiop,' \n');
fprintf(fiop,formatesp, size(A),mean(A),std(A));

formatesp3=['Number of useful points 1-2, average cos, stdesv, %d, %d, %5.3f, %5.3f \n'];
formatesp4=['promecos_2, averageradcos, averageradcos2, averageradcossin 1-2 %5.3f, %5.3f, %5.3f, %5.3f, %5.3f \n'];

fprintf(fiop,' \n');

fprintf(fiop,'Data cos-weighted and rad-weighted A1 \n');
fprintf(fiop,' \n');
%first in case of 'N', here only one file
[numRows,numCols] = size(A1);
m=numRows;
n=numCols;
%print in screen to see it works
disp('m, n from A1: ');
disp(m);
disp(n);
np=0;
sumecsecs=0;
sumdes=0;
e2=0;
e1=0;
ee=0;
eea=0;
eea2=0;
eea2s=0;
eeasis=0;
sumecsecs2=0;
sumecsecs2s=0;
mm=1;

%Attention separated by hemispheres

for k=1:1:floor(m/2)   
    ecsecs=abs(sin(floor(k)*pi/m));
    ecsecs2=ecsecs*ecsecs;
    ecsecssi=sqrt((1.-ecsecs2));
    for l=1:1:n
        %FER AMB if() end
        if (~isnan(A1(k,l)))
            % if A1(k,l) different from nan
            A11(mm)= A1(k,l);
            A111(mm)=ecsecs;
            mm=mm+1;
            e0=ecsecs*A1(k,l); %temperature in ºC
            e1=e1+e0;
            e2=e2+e0*A1(k,l);
            eea= (A1(k,l)+273.15)^4; %radiation
            eea2s=eea2s+eea*ecsecs;
            eeasis=eeasis+eea*ecsecs*ecsecssi;
            ee= ee+eea*ecsecs2;
            sumdes= sumdes+eea*eea*ecsecs; % normalizations
            sumecsecs= sumecsecs+ecsecs;
            sumecsecs2=sumecsecs2+ecsecs2;
            sumecsecs2s=sumecsecs2s+ecsecs*ecsecssi;
            np=np+1;
        else
            ;
        end;
    end;
end;

sumecsecs2s2=0;  %Careful!
eeasis2=0;
np2=0;
for k=floor(m/2)+1:1:m   % second hemisphere
%    ecsecs=cos((1-2*(k-.5)/m)*pid2);
    ecsecs=abs(sin(floor(k)*pi/m));
    ecsecs2=ecsecs*ecsecs;
    ecsecssi=sqrt((1.-ecsecs2));
    for l=1:1:n
        % if() - end
        if (~isnan(A1(k,l)))
            % if A1(k,l) different from nan %temperature ºC
            A11(mm)= A1(k,l);
            A111(mm)=ecsecs;
            mm=mm+1;
            e0=ecsecs*A1(k,l);
            e1=e1+e0;
            e2=e2+e0*A1(k,l); %radiation
            eea= (A1(k,l)+273.15)^4;
            eea2s=eea2s+eea*ecsecs;
            eeasis2=eeasis2+eea*ecsecs*ecsecssi;
            ee= ee+eea*ecsecs2;
            sumecsecs= sumecsecs+ecsecs;
            sumecsecs2=sumecsecs2+ecsecs2;
            sumecsecs2s2=sumecsecs2s2+ecsecs*ecsecssi;
            np2=np2+1;
           else
           
        end;
    end;
end;

promecos=e1/sumecsecs;
desvest2=sqrt(e2/sumecsecs-promecos*promecos);%sqrt
% Attention come back from radiation to temperature
promecos_2=sqrt(sqrt(ee/sumecsecs2));
primerm=promecos;
primerdesv=desvest2;
%sumdes2=sqrt(sqrt(sqrt((sumdes/sumecsecs2-ee*ee/sumecsecs/sumecsecs))));
disp('Average taking into account cos(theta)');
disp('promecos: ');
disp(promecos);
%disp('desvestcos: ');
%disp(sumdes2);
%fprintf(fiop,formatesp3,np,np2,promecos,desvest2,promecos_2,sumdes2,(ee/sumecsecs)^.25,(eea2s/sumecsecs2)^0.25,(eeasis/sumecsecs2s)^.25,(eeasis2/sumecsecs2s2)^.25);
fprintf(fiop,formatesp3,np,np2,promecos,desvest2);
fprintf(fiop,formatesp4, promecos_2-273.15, (ee/sumecsecs)^.25-273.15,(eea2s/sumecsecs)^0.25-273.15, ...
    (eeasis/sumecsecs2s)^.25-273.15,(eeasis2/sumecsecs2s2)^.25-273.15);
;
fprintf(fiop,' \n');

%computation spatial contrast
sc0=0.0;
sc00=0.0;
sumecs=0.0;
mp=0;
for k=2:1:m-1 %files
    %weight depends on file number! multiply sum by the cosinus
      ecsecs=abs(sin(floor(k)*pi/m));
      sc=0.0;
np=0;
%Do only if no NAN!!!

    for l=2:1:n-1 %each column    compute Sum[abs(value-value-1)+...                
if (~isnan(A1(k,l)) && ~isnan(A1(k,l-1)) && ~isnan(A1(k,l+1)) && ~isnan(A1(k-1,l)) && ~isnan(A1(k+1,l)))
  sc=sc+(abs(A1(k,l)-A1(k,l-1))+abs(A1(k,l)-A1(k,l+1))+abs(A1(k,l)-A1(k-1,l))+abs(A1(k,l)-A1(k+1,l)))/4;
  np=np+1;
  np;
  %
else
    %don't do anything
    ;
end; %end of if
;
    end; %end of loop l (n)
;
if (np>0)
 sc=sc*ecsecs/np %weight each file and normalize
 sumecs=sumecs+ecsecs;
 sc0=sc0+sc;   %Cumulative for files danger if it's sc=0 don't count
else
    ;
end %of np>0
;
end; %of loop k (m)
sc0=sc0/sumecs;
sc00=sc0/(maxI1-minI1);
; %divide by (m-2)->sumcos and by abs(maxI1-minI1) %danger assumption 
;%print the result
fprintf(fiop,' \n');
fprintf(fiop,'Spatial contrast (in units of deltacontiguous divided by number, and by Tmax-Tmin) %f %f \n',sc0,sc00);
fprintf(fiop,' \n');
;
%txt=input('Enter text: ','s');
%stop if needed
;

% define vector weight (sinus) and do the reshape equal, then output
% the average and standard deviations attention, to compute with cos is needed
%to go to matrix
disp('A: elements, average, std: ');
disp(size(A));
disp(mean(A));
disp(std(A));
;
%Attention

figure(2); %figure of violins
nexttile;
Z=[A];
group1=[ones(size(A))];
v=violinplot(group1,Z);
clear ylim;
ylim([minI1 maxI1]);
ylabel("Lineal distribution of Radiant Temperature in ºC");
% xticklabels({' ','7 h', ' ', '11 h', ' ', '15 h',' ','19 h'});
xticklabels({' ',' ',' ',' ',' ','15 h'});
xlabel("Revolucion 07072026");
;
% annotate sample sizes /Better mean and stdev
%ns = [numel(A), numel(B), numel(Bu), numel(C), numel(Cu), numel(Cv), numel(D)];
%as = [primerm, segonm, tercerm, quartm];
%ds = [primerdesv, segondesv, tercerdesv, quartdesv];
as = [primerm];
ds = [primerdesv];
for k=1:1 % before it was 1... 4
    text(k,min(A)+9,sprintf('av= %.2f',as(k)),... 
        'HorizontalAlignment','center');
    text(k,min(A)+2,sprintf('std= %.2f',ds(k)),...
        'HorizontalAlignment','center');
end;
;
%Violinplot cos corrected
%nexttile;
%Z=[A11(:); B11(:); C11(:); D11(:)];
%ZZ=[A111(:); B111(:); C111(:); D111(:)];
%group1=[ones(numel(A11),1); ones(numel(B11),1)+1; ones(numel(C11),1)+2; ones(numel(D11),1)+3];

clear ylim;
figure(3);
DensityScale='area';
minI1=-35;
maxI1=85;
[densVal1, evalPts1]=ksdensity(A11, 'Weights', A111);
v1=violinplot(EvaluationPoints=evalPts1, DensityValues=densVal1);
ylabel("Distribution (cos-w) of Radiant Temperature in ºC");
ylim([minI1 maxI1]);
d1max=max(densVal1);
xticklabels({sprintf('Width of distribution %.5f',2*d1max)});
xlabel("Revolucion 15 h 07072026");
;
% Decomposition in spherical harmonics
% function a1 = l1SphHarmCoeffs(A1)
[m, n] = size(A1);

% Grid on the sphere
theta = ((1:m)' - 0.5) * pi / m;      % 0..pi
phi   = ((0:n-1) + 0.5) * 2*pi / n;    % 0..2pi

[Phi, Theta] = meshgrid(phi, theta);
dOmega = sin(Theta) * (pi / m) * (2*pi / n);

% Remove NaNs
T = A1;
valid = ~isnan(T);
T(~valid) = 0;

% l = 1 spherical harmonics
Y1m1 =  sqrt(3/(8*pi)) * sin(Theta) .* exp(-1i * Phi);
Y10  =  sqrt(3/(4*pi))  * cos(Theta);
Y11  = -sqrt(3/(8*pi)) * sin(Theta) .* exp( 1i * Phi);

% Coefficients a_lm = integral T * conj(Y_lm) dOmega
a1m1 = sum(T(valid) .* conj(Y1m1(valid)) .* dOmega(valid));
a10  = sum(T(valid) .* conj(Y10(valid))  .* dOmega(valid));
a1mas1  = sum(T(valid) .* conj(Y11(valid))  .* dOmega(valid));

%a1 = struct('m_1', a1m1, 'm0', a10, 'p1', a1mas1);
fprintf(fiop,' \n');
fprintf(fiop,'Case 15 h Revolucion 07072026 size m*n %.5d %.5d \n',m,n);
fprintf(fiop,' \n');
fprintf(fiop,'Spherical Harmonics Ylm (complex values) for l=1, m=-1, 0, 1 : %f %f %f %f %f %f \n',real(a1m1), imag(a1m1),real(a10),imag(a10),real(a1mas1),imag(a1mas1));
fprintf(fiop,' \n');
% coefficients:
% a1m1 = Y(1,-1)
% a10  = Y(1, 0)
% a1mas1  = Y(1, 1)

v = [a1m1 - a1mas1;
    -1i*(a1m1 + a1mas1);
    sqrt(2)*a10];

% for a real field, this should be (almost) real
v = real(v);
v = v / norm(v);

theta = acos(max(-1,min(1,v(3))));   % 0..pi
phi   = atan2(v(2), v(1));           % -pi..pi
phi   = mod(phi, 2*pi);              % reduce to first circ.
fprintf(fiop,' \n');
fprintf(fiop,'Preferred direction of temperature according to first Spherical Harmonics \n');
fprintf(fiop,' \n');
fprintf(fiop,'Direction in xyz (vector), theta, fi (rad): %f %f %f  %f %f \n',v(1), v(2), v(3), theta, phi);
fprintf(fiop,' \n');
fprintf(fiop,'Direction in theta, fi in º: %f %f \n', theta*180/3.1415926, phi*180/3.1415926);
;

