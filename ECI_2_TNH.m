%% Coordinate frame transformation: ECI -> TNH
function R = ECI_2_TNH(r,v)
    ix = v/norm(v);
    iz = cross(r,v);
    iz = iz/norm(iz);
    iy = cross(iz,ix);
    R = inv([ix,iy,iz]);
end