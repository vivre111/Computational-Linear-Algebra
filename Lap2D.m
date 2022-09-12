function [A, b] = Lap2D(m)
    e = ones(m^2,1);
    A = spdiags([-e, -e, 4*e, -e, -e],[-m,-1,0,1,m], m^2, m^2);
    for i = 1:(m-1)
        for j = 1:(m-1)
            A(m*i, m*j+1) = 0;
            A(m*i+1, m*j) = 0;
        end
    end
    h = 1/ (m+1);
    b = zeros(m^2,1);
    x_low = ceil(0.25/h);
    x_high = floor( 0.45/h );
    y_low = ceil(0.5/h);
    y_high = floor(0.7/h);
    for i = x_low:x_high
        for j = y_low:y_high
            b((j-1)*m+i)=1;
        end
    end
    x_low = ceil(0.7/h);
    x_high = floor( 0.9/h );
    y_low = ceil(0.15/h);
    y_high = floor(0.35/h);
    for i = x_low:x_high
        for j = y_low:y_high
            b((j-1)*m+i)=1;
        end
    end
end

