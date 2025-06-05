function y = prbs(n,m)

% PRBS  Pseudo-Rausch-Binärsignal
% Es wird ein PRBS mit den Werten +1 und -1 erzeugt; die Länge
% n des Schieberegisters kann zwischen 1 und 11 liegen. Jeder
% Signalwert wird m Schritte gehalten. Daraus resultiert eine
% Signallänge von m*(2^n-1).
%
% Syntax: y = prbs(n);
%         y = prbs(n,m);
% Input : n : Anzahl Registerstufen
%         m : Zeiteinheiten pro Takt (Default: m=1)
% Output: y : PRBS (Spaltenvektor), Länge m*(2^n-1)

% LOLIMOT Nonlinear System Identification Toolbox, Version 2.1
% Oliver Nelles and Alexander Fink, 13-May-2002
% Institute of Automatic Control, Darmstadt University of Technology, Germany
% Copyright (c) 2002 by Oliver Nelles and Alexander Fink


if nargin==1, m=1; end                 % Default-Wert für m
z=-ones(n,1); 
y=zeros(m*(2^n-1),1);    % Initialisierung
s={[1 1],[1 2],[2 3],[3 4],[3 5], ...  % Registerstufen
   [5 6],[4 7],[1 2 7 8],[5 9],[7 10],[9 11]};
for k=1:m:length(y)
  z=[prod(z(s{n})); z(1:n-1)];         % Register schieben
  y(k:k+m-1)=z(n)*ones(m,1);           % PRBS auslesen
end
