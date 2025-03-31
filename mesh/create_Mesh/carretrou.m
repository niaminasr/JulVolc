function [x,y] = carretrou(bs,s)
% Create a unit circle centered at (0,0) using four segments.
switch nargin
    case 0
        x = 8; % four edge segments
        return
    case 1
        A = [0, 0, 0, 0, 0, pi/2, pi, 1.5*pi ;    % start parameter values
             1, 1, 1, 1, pi/2, pi, 1.5*pi, 2*pi; % end parameter values
             1, 1, 1, 1, 0, 0, 0, 0 ;    % region label to left
             0, 0, 0, 0, 1, 1, 1, 1];    % region label to right
        x = A(:,bs); % return requested columns
        return
    case 2
        taillebs = size(bs);
        tailles = size(s);

        % Parameters of the circle

        xc = 0.5;
        yx = 0.5;
        radius = 0.1;

        if (taillebs == [1 1]) % Ici un seul segment frontière et plusieurs s
            switch(bs)
                case 1
                    x = s;
                    y = zeros(size(s));
                case 2
                    x = ones(size(s));
                    y = s;
                case 3
                    x = 1-s;
                    y = ones(size(s));
                case 4
                    x = zeros(size(s));
                    y = 1-s;
                case 5
                    x = xc + radius*cos(s);
                    y = xc + radius*sin(s);
                case 6
                    x = xc + radius*cos(s);
                    y = xc + radius*sin(s);
                case 7
                    x = xc + radius*cos(s);
                    y = xc + radius*sin(s);
                case 8
                    x = xc + radius*cos(s);
                    y = xc + radius*sin(s);
            end
           
            
        else % Ici s et bs ont la même taille
            x = zeros(size(s));
            y = x; % Initialize with zeros
            for index = 1:tailles(1)
                for jndex = 1:tailles(2)
                    switch(bs(index,jndex))
                        case 1
                            x(index,jndex) = s(index,jndex);
                            y(index,jndex) = 0;
                        case 2
                            x(index,jndex) = 1;
                            y(index,jndex) = s(index,jndex);
                        case 3
                            x(index,jndex) = 1-s(index,jndex);
                            y(index,jndex) = 1;
                        case 4
                            x(index,jndex) = 0;
                            y(index,jndex) = 1-s(index,jndex);
                        case 5
                            x(index,jndex) = xc + radius*cos(s(index,jndex));
                            y(index,jndex) =xc + radius*sin(s(index,jndex));
                        case 6
                            x(index,jndex) = xc + radius*cos(s(index,jndex));
                            y(index,jndex) =xc + radius*sin(s(index,jndex));
                        case 7
                            x(index,jndex) = xc + radius*cos(s(index,jndex));
                            y(index,jndex) =xc + radius*sin(s(index,jndex));
                        case 8
                            x(index,jndex) = xc + radius*cos(s(index,jndex));
                            y(index,jndex) =xc + radius*sin(s(index,jndex));
                    end
                end
            end
        end
end