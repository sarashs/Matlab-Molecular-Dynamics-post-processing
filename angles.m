function [ angle ] = angles( atom1atom2,atom1atom3)
%Calculates the atom2-atom1-atom3 angle
% atomi is a [x y z] array showing the location of an atom
a=atom1atom2;
b=atom1atom3;
angle=180*acos(a*b'/(norm(a)*norm(b)))/pi;
end

