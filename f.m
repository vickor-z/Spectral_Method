function [ fu ] = f( x )
fu = (12.*cos(x).*(5-4.*cos(x))-96.*sin(x).^2)./((5-4.*cos(x)).^3);
end

