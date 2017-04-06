# geometry-processing-project
Conformal Quaternion Based Mesh (Re)Construction is a method for deforming a seed mesh into the target mesh implied by a 3D point cloud. The images below represent a sphere transformed into a human head. For more details see /abstract/technical-abstract.pdf
<table>
  <thead>
  </thead>
  <tbody>
    <tr>
      <td>
        <img src="/presentation/images/front_sphere_to_head_overnight.PNG" alt="" width="300px" height="450px">
      </td>
      <td>
        <img src="/presentation/images/side_sphere_to_head_overnight.PNG" alt="" width="300px" height="450px">
      </td>
      <td>
        <img src="/presentation/images/back_sphere_to_head_overnight.PNG" alt="" width="300px" height="450px">
      </td>
    </tr>
  </tbody>
</table>

Setup

- Clone the repo to your local machine
- Create a build folder alongside main.cpp
- Run cmake from one level above the build folder with the target being the build folder. Be sure to use a 64 bit compiler as the large matrices involved in this project cannot be handled by 32 bit. You will be runtime memory errors otherwise. 
