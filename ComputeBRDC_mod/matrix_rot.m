function R=matrix_rot(teta,verso,asse)
%
%CHIAMATA: R=matrix_rot(teta,verso,asse)
%la function matrix_rot costruisce matrici di rotazione
%in input angolo di rotazione, verso e asse di rotazione
%
%utilizzo: 
%inserire teta in radianti
%
%verso='O' per rotazioni orarie,
%verso='A' per rotazioni antiorarie
%
%asse='x' per rotazioni intorno all'asse x
%asse='y' per rotazioni intorno all'asse y
%asse='z' per rotazioni intorno all'asse z 

if verso=='A'
    teta=-teta;
elseif verso=='O'
else
    error('errore inserimento verso rotazione')
end

if asse=='x'
   R=[1 0 0;0 cos(teta) -sin(teta);0 sin(teta) cos(teta)];
elseif asse=='y'
   R=[cos(teta) 0 sin(teta);0 1 0;-sin(teta) 0 cos(teta)];
elseif asse=='z'
   R=[cos(teta) -sin(teta) 0;sin(teta) cos(teta) 0;0 0 1];
else
   error('errore inserimento asse rotazione') 
end