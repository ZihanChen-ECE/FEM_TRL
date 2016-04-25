clc
clear
gmsh_filename = 'rect.msh';
fid = fopen(gmsh_filename);

 while ( 1 )

    text = fgetl ( fid )
    class(text)
    length(text)
    if ( text(end) == -1 )
      break
    end
 end