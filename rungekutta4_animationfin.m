clear all;
close all;
clc;

function [resultat]=somme(W,X,h)
  %W positions des aimants 1 à N
  N=length(W); %nombre d'aimants
  resultat=zeros(2,1);
  for i=1:N
    resultat(1)=resultat(1) + (W(i,1)-X(2,1))/(( ((W(i,1)-X(2,1)).^2 +(W(i,2)-X(2,2)).^2).^2  .+(h^2)).^(5/2));
    resultat(2)=resultat(2) + (W(i,2)-X(2,2))/(( ((W(i,1)-X(2,1)).^2 +(W(i,2)-X(2,2)).^2).^2  .+(h^2)).^(5/2));
  endfor
end

function [A] = f(X,W,t,h,b)
  A = zeros(2,2);
  S=somme(W,X,h);
  A(1,:)=-b.*X(1,:) - X(2,:) + S';
  A(2,:)=X(1,:);
end

function [U,V]= RK(pas,T,X0,h,b,W)
  %U : vecteur position du pendule, V vecteur vitesse du pendule
  X=X0; %X0 CI matrice 2*2
  n=T/pas; %nb d'itérations
  U=[X0(2,:)'];
  V=[X0(1,:)'];
  for i=1:n
    k1=f(X,W,i*pas,h,b);
    k2=f(X+pas*k1/2,W,(i+1)*pas/2,h,b);
    k3=f(X+pas*k2/2,W,(i+1)*pas/2,h,b);
    k4=f(X+pas*k3,W,(i+1)*pas,h,b);
    X=X+ (pas/6)*(k1+2*k2+2*k3+k4);
    U=[U X(2,:)']; %positions
    V=[V X(1,:)']; %vitesses
  endfor
  
  end




h = 0.5 ;
b = 0.1;
N = 3; % Nombre d'aimants
d = 1; %distance des aimants à l'origine
pas = 0.01; %pas temporel
T   = 60;  %Durée de l'expérience 
X0  = [ 0 0; 0 1.01]; %Condition initiale matrice de taille 2x2, [x' y'; x y] 
W   = [];
hold on;
for i=1:N
  xi=cos(i*2*pi/N)*d;
  yi=sin(i*2*pi/N)*d;
  X = [xi yi];
  W = [ W; X];
  plot(xi,yi,'.','Markersize',20,'Color','b'); %on trace les aimants
 end




%plot(U(1,:),U(2,:)) ;
%X2=[0 0; 2*cos(2*pi/N)*d 2*sin(2*pi/N)*d+0.001];


X2  = [ 0 0; 2 2];
[U2,V2] = RK(pas,T,X2,h,b,W);
%plot(U2(1,:),U2(2,:)) ;

for k=1:length(U2(1,:))
  plot(U2(1,k),U2(2,k));
 pause(.0000000000000001)
end

 
figure(2);
temps=[];
for i=1:length(U)
  temps=[temps pas*i];
  endfor

%plot(temps,U(1,:));