%
%                     SIMMULACIÓN DE UN SISTEMA PLANETARIO
%
% Santiago Fernández Gutiérrz Zamora                                 
%==========================================================================

% PARÁMETROS

%        Datos de wikipedia
% ( Nombre , masa [kg] , radio [km] , radio orbital [km] , angulo inclinación [deg], velocidad orbital [km/s])   
datos = {'Sol' , 1.989e30 , 695508 , 0 , 0 , 0 ;...
         'Mercurio' , 3.302e23 , 2439.7 , 57894376 , 7 , 47.8725 ;...
         'Venus' , 4.869e24 , 6051.8 , 108208930 , 3 , 35.0214 ;...
         'Tierra' , 5.9736e24 , 6371 , 149597870.691 , 7 , 29.78 ;...
         'Luna' , 7.349e22 , 1737.1 , 149597870.691+384403 , 7 , 29.78+1.019927755507821 ;...
         'Saturno V' , 2.8e6 , 2e-3 , 4.224400625146314e+04+149597870.691 , 7 , 3.072070826337914+29.78 ;... % la nave
         'Marte' , 6.4171e23 , 3389.5 , 249200000 , 5 , 24.007 }; %...

% Indice del movil que se fija al centro [1 es el sol]
movil_fijo = 1;

% parámetros de la nave
hay_galletazo = 1; % [1 : true] [0 : false]
nav_indx = 6; % indice de la nave con respecto a los otros cuerpos
t_galletazo = 4000; % el subintervalo en el cual se aplica el galletazo 0 < t_g < Nk
pcm = 0.7; % porcentaje de la masa que se consume 0 < c < 1
ve = 2.4e3; % exhaust velocity (wikipedia para saturno V)
u_galletazo = [-1/sqrt(3);1/sqrt(3);-1/sqrt(3)]; % vector de la dirección del galletazo


% subintervalo de tiempo (en segundos)
deltaT = 3600;
% 60 : un minuto || 3600 : una hora || 86400 : un día

% tiempo total (en segundos) tiene que ser un múltiplo del subintervalo
tt = deltaT*8760*2;
% tt = 28*24*60*60; % un ciclo lunar -ish

Nk = tt/deltaT; % número de subintervalos 

G = 6.67408E-11; %  constante gravitacional :: m3 kg-1 s-2  

% INSTANCIAR LOS CUERPOS

[nc,ln] = size(datos);
cuerpos = Movil.empty(nc,0);
r = zeros(3,nc);
v = zeros(3,nc);
pos = zeros(3,nc,Nk);  
for k = 1:nc
    cuerpos(k).nombre = datos(k,1);
    cuerpos(k).masa = cell2mat(datos(k,2));
    cuerpos(k).radio = cell2mat(datos(k,3))*1000; % a metros
    cuerpos(k).dist_sol = cell2mat(datos(k,4))*1000;
    cuerpos(k).ang_inc = cell2mat(datos(k,5))*pi/180; % radianes
    cuerpos(k).vel_orb = cell2mat(datos(k,6))*1000; % m/s
    cuerpos(k) = cuerpos(k).trayectoria(Nk);
    r(:,k) = cuerpos(k).pos_inicial';
    pos(:,k,1) = r(:,k);
    v(:,k) = cuerpos(k).vel_inicial';
end

% ITERAR LA SOLUCIÓN

for t = 1:Nk 
    
    % CALCULO DE FUERZAS
    
    dr= zeros(nc,nc,3); % vectores de diferencias
    
    dr2 = zeros(nc); % distancias al cuadrado
    
    u  = zeros(nc,nc,3); % Vectores unitarios de dirección de fuerzas
    
    F = zeros(nc,nc,3); % Vectores de fuerza entre cada movil
    
    Fres = zeros(3,nc); % vector de fuerza resultante (por cada movil)
    
    a = zeros(3,nc); % vector de aceleración (por cada movil)
    
    for j = 1:nc % columnas
        for i = 1:nc % renglones
            if j ~= i
                dr(i,j,:) = r(:,i) - r(:,j);
                if j<i
                    dr2(i,j) = sum(dr(i,j,:).^2); 
                    u(i,j,:) = dr(i,j,:)/sqrt(dr2(i,j)); 
                    F(i,j,:) = ((G*cuerpos(j).masa*cuerpos(i).masa)/dr2(i,j))*u(i,j,:); 
                else
                    dr2(i,j) = dr2(j,i);
                    u(i,j,:) = -u(j,i,:);
                    F(i,j,:) = -F(j,i,:);
                end
            end
        end
        
        Fres(:,j) = sum(F(:,j,:)); % Fuerza resultante
        
        a(:,j) = Fres(:,j)/cuerpos(j).masa; % Aceleración
        
    end
    
    % EULER
    
    v = v + a * deltaT; % Velocidad
    
    if hay_galletazo == 1 && t == t_galletazo
       dm = cuerpos(nav_indx).masa - cuerpos(nav_indx).masa*pcm; % dry mass
       [dv,cuerpos(nav_indx)] = cuerpos(nav_indx).galletazo(ve,dm);
       v(:,nav_indx) = v(:,nav_indx) + dv*u_galletazo;
    end
  
    r = r + v * deltaT; % Posición
    
    for k = 1:nc
        pos(:,k,t+1) = r(:,k);
    end
    
end

pos2 = zeros(3,nc,Nk+1);
for ca = 1:nc
    if ca ~= movil_fijo
        for coo = 1:3
            pos2(coo,ca,:) =  pos(coo,ca,:) - pos(coo,movil_fijo,:);
        end
    end
end

% GUARDAR LAS TRAYECTORIAS EN CADA CUERPO
for k = 1:nc
   cuerpos(k).tryct(:,:) = reshape(pos2(:,k,:),3,[])'; 
end



% GRÁFICAS

numlgnd = nc;
if hay_galletazo == 1
   pos_galletazo = cuerpos(nav_indx).tryct(t_galletazo,:); % la posición donde se acelera
   numlgnd = numlgnd + 1; 
end
legends = string.empty(numlgnd,0);
plot3(cuerpos(1).tryct(:,1),cuerpos(1).tryct(:,2),cuerpos(1).tryct(:,3),'o') % el sol
legends(1) = cuerpos(1).nombre;
hold on
grid on
for k = 2:nc
   %     ( x , y , z )
   plot3(cuerpos(k).tryct(:,1),cuerpos(k).tryct(:,2),cuerpos(k).tryct(:,3),'.')
   legends(k) = cuerpos(k).nombre;
end
if hay_galletazo == 1
    plot3(pos_galletazo(1),pos_galletazo(2),pos_galletazo(3),'+');
    legends(k+1) = 'Galletazo';
end
title('Sistema solar');
xlabel('Coordenada X (mts)');
ylabel('Coordenada Y (mts)');
zlabel('Coordenada Z (mts)');
legend(legends);
hold off