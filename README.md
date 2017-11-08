# pathtracer
Minimalist Quick'n Dirty Path Tracer. No keyboard needed, just drag'n drop files.

## Handles:
* Phong BRDF model ; some defaults materials from Ngan's fits
* Environment maps
* Normal maps, Alpha maps, textures for Albedo, Specular Ks, Shininess, Refraction, Refraction Index
* OBJ (with materials) and VRML import, XYZ import for point clouds (can also estimate normals)
* Display of edges, supports non-triangular faces
* Topological information
* Filtering, DoF, one spherical light source, Fresnel reflection on transparent objects.
* Limited and buggy support of fog.

## Images:

<table>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/bot.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/botNew.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/botTransp.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 2.5 million triangles</td>
  <td> 2.5 million triangles</td>
  <td> 2.5 million triangles, entirely transparent</td>
  </tr>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/babyroom.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/americanGirl.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/theBuilding.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 1.8 million triangles</td>
  <td> 1.6 million triangles</td>
  <td> 3.1 million triangles, includes transparency</td>
  </tr>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/ponies.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/man.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/lion.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 133 k triangles</td>
  <td> 260 k polygons</td>
  <td> 1.8 k triangles + normal map</td>
  </tr>
   <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/antiqueOffice.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/antiqueRoom.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/ship.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 23.7 million triangles</td>
  <td> 1.8 million triangles</td>
  <td> 2.5 k polygons + normal map, alpha map and edges display</td>
  </tr>
  
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/alien.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/mepointcloud.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 240 k triangles</td>
  <td> 4.5 millions points Point Cloud</td>
  </tr>
  </table>
