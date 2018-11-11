function [ x, y, z ] = subtract( atom1, atom2, maxx, maxy, maxz )
%atomi= [x,y,z]
% gives a vector atom1-atom2
       if abs(atom1(1)-atom2(1))<(maxx/2)
           x=atom1(1)-atom2(1);
       elseif atom1(1)-atom2(1)<0
           x=maxx+atom1(1)-atom2(1);
       else
           x=-maxx+atom1(1)-atom2(1);
       end
      
       if abs(atom1(2)-atom2(2))<(maxy/2)
           y=atom1(2)-atom2(2);
       elseif atom1(2)-atom2(2)<0
           y=maxy+atom1(2)-atom2(2);
       else
           y=-maxy+atom1(2)-atom2(2);
       end

       if abs(atom1(3)-atom2(3))<(maxz/2)
           z=atom1(3)-atom2(3);
       elseif atom1(3)-atom2(3)<0
           z=maxz+atom1(3)-atom2(3);
       else
           z=-maxz+atom1(3)-atom2(3);
       end       
end

