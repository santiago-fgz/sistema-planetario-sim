classdef Movil
    % clase Movil
    
    % nombre: string
    % masa: [kg]
    % radio: radio del planeta [m]
    % dist_sol: distancia media del sol [m]
    % ang_inc: angulo de inclinación [deg]
    % vel_orb: velocidad orbital media [km/s]
    
    % Nk: número de subintervalos
    % tryct: [ xi , yi , zi ] [m] :: matriz con la trayectoria
    
    
    properties
        nombre
        masa {mustBeNumeric}
        radio {mustBeNumeric}
        dist_sol {mustBeNumeric}
        ang_inc {mustBeNumeric}
        vel_orb {mustBeNumeric}
        Nk {mustBeNumeric}
        tryct {mustBeNumeric}
    end
    
    methods
        
        function obj = movil(name,mass,radius,sun_dist,inc_ang,orb_vel)%,strtPos,strtVel)
            % CONSTRUCTOR
            if nargin == 4 
                obj.nombre = name;
                obj.masa = mass;
                obj.radio = radius;
                obj.dist_sol = sun_dist;
                obj.ang_inc = inc_ang*pi/180;
                obj.vel_orb = orb_vel;
            end
        end
        
        function obj = trayectoria(obj,n)
           % definir el número de intervalos para la trayectoria
           % m = m.trayectoria(1000);
           if nargin == 2
              obj.Nk = n+1;
              obj.tryct = zeros(obj.Nk,3); 
              %obj.tryct(1,:) = obj.pos0; 
           end
        end
        
        function pos = pos_inicial(obj)
            xPos = obj.dist_sol*cos(obj.ang_inc);
            yPos = 0;
            zPos = obj.dist_sol*sin(obj.ang_inc);
            pos = [xPos yPos zPos];
        end
        
        function vel = vel_inicial(obj)
            vel = [0 obj.vel_orb 0];
        end
        
        function [dv,obj] = galletazo(obj,ve,dm)
            dv = ve*log(obj.masa/dm);
            obj.masa = dm;
        end
        
    end
end