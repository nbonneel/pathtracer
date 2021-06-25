# pathtracer
Minimalist Quick'n Dirty Path Tracer. (almost) No keyboard needed, just drag'n drop files.

## Supports:
* Phong BRDF model ; some defaults materials from Ngan's fits, as well as support for measured BRDFs from Matusik
* Environment maps
* Normal maps, Alpha maps, textures for Albedo, Specular Ks, Shininess, Refraction, Refraction Index
* OBJ (with materials) and VRML import, XYZ import for point clouds (can also estimate normals)
* Display of edges, supports non-triangular faces
* Topological information
* Filtering, DoF, one spherical light source, Fresnel reflection on transparent objects.
* Rendering for lenticular images (can buy sheets at vuethru.com or dplenticular.com) and camera arrays (launch computation with the Animation offline render button -- sorry for that ; it will directly write all your files on disk).
* Limited and buggy support of fog.
* Some support for integrating photos (background photo, "ghost" objects that receive shadows and indirect light but are not displayed)

## Compilation:
* Dependency: Only wxWidget 3.1.0 (to download), CImg (included) and Nanoflann (included). Optional: Embree. Deactivate by removing #define USE_EMBREE. If using Embree, compile Embree with EMBREE_RAY_MASK activated (deactivated by default). Embree roughly halves the intersection time.
* Windows: use .sln solution for Visual Studio. Edit path to wxWidget.
* Linux/MacOS: use CMake, should work too.

## Usage:
* Drag'n drop OBJ, VRML, XYZ (pointsets) files onto window. Also supports .Yarn files (generated by [https://www.cs.cornell.edu/projects/ctcloth/#proc-sig15](this paper)), Matusik's isotropic .binary BRDF files ([https://cdfg.csail.mit.edu/wojciech/brdfdatabase](BRDF dataset here)), .Titopo BRDFs ([https://perso.liris.cnrs.fr/guillaume.lavoue/data/BRDFs/](from here)), .seg segmentation files ([https://segeval.cs.princeton.edu/](from here, for example) -- one segment ID per face) and .lab segmentation files (where did I get these from?!).
* Mouse down (left, middle, right) for camera control, mouse wheel for moving forward
* Click to select object ; see the cursor depth in the status bar (useful to change the focus distance)
* Shift+Mouse to move/rotate currently selected object
* Drag'n drop texture files or envmap files onto their slot on the right panel, or use the popup menu to change their values.
* If materials are changed, they can be re-exported into a new .MTL file (menu File)
* Menu Info gives topological info on the currently selected mesh (only works for meshes, not point clouds nor the default parametric geometries)
* Scene can be saved in an .scn file

## Images:

<table>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/bot.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/botNew.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/botTransp.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 2.5 million triangles <a href="https://www.cgtrader.com/free-3d-models/character/sci-fi/24h-burnoutbot">(src)</a></td>
  <td> 2.5 million triangles <a href="https://www.cgtrader.com/free-3d-models/character/sci-fi/24h-burnoutbot">(src)</a></td>
  <td> 2.5 million triangles, entirely transparent <a href="https://www.cgtrader.com/free-3d-models/character/sci-fi/24h-burnoutbot">(src)</a></td>
  </tr>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/babyroom.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/americanGirl.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/theBuilding.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 1.8 million triangles <a href="https://www.cgtrader.com/free-3d-models/interior/bedroom/decor-for-children">(src)</a></td>
  <td> 1.6 million triangles <a href="https://www.turbosquid.com/3d-models/free-obj-model-american-beauty/904057">(src)</a></td>
  <td> 3.1 million triangles, includes transparency <a href="https://www.blendswap.com/blends/view/73806">(src)</a></td>
  </tr>
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/ponies.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/man.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/lion.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 133 k triangles <a href="https://www.blendswap.com/blends/view/70960">(src)</a></td>
  <td> 260 k polygons <a href="https://www.cgtrader.com/free-3d-models/character/man/3d-scan-man">(src)</a> </td>
  <td> 1.8 k triangles + normal map <a href="https://www.cgtrader.com/free-3d-models/animals/mammal/lion-statue-low-poly">(src)</a></td>
  </tr>
   <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/antiqueOffice.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/antiqueRoom.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/ship.jpg" width="300"> </td>
  </tr>
  <tr> 
  <td> 23.7 million triangles <a href="https://www.blendswap.com/blends/view/83895">(src)</a></td>
  <td> 1.8 million triangles <a href="https://www.blendswap.com/blends/view/86774">(src)</a></td>
  <td> 2.5 k polygons + normal map, alpha map and edges display <a href="https://www.cgtrader.com/free-3d-models/watercraft/other/ghost-ship">(src)</a></td>
  </tr>
  
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/alien.jpg" width="300"> </td>
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/mepointcloud.jpg" width="300"> </td>
  <td><a href="https://www.youtube.com/watch?v=_L8Au4MOjr8"><img src="https://img.youtube.com/vi/_L8Au4MOjr8/0.jpg" width="300"></a></td>  
  </tr>
  <tr> 
  <td> 240 k triangles <a href="https://www.cgtrader.com/free-3d-print-models/miniatures/figurines/raisher-the-sky-reaper-printable-bust">(src)</a></td>
  <td> 4.5 millions points Point Cloud</td>
  <td> Lenticular images</td>
  </tr>
  
  <tr>  
  <td> <img src="https://github.com/nbonneel/pathtracer/raw/master/img/lionPhoto.jpg" width="300"> </td>
  <td> </td>
  <td> </td>  
  </tr>
  <tr> 
  <td>Lion on a ghost plane + Background photo</td>
  <td> </td>
  <td> </td>
  </tr>  
  
  </table>
