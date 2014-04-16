HOMEWORK 3: RAY TRACING, RADIOSITY, & PHOTON MAPPING

NAME:  Andrew Leing

Compiled on Windows 7 using GLFW
Ran using Visual Studio Express 2013

TOTAL TIME SPENT:  38 hours
Please estimate the number of hours you spent on this assignment.


COLLABORATORS: 
You must do this assignment on your own, as described in the 
Academic Integrity Policy.  If you did discuss the problem or errors
messages, etc. with anyone, please list their names here.

-Got assistance from office hours

RAYTRACING:
For the first reflective spheres, shadows are generated for each one even though the homework page doesn't have them enabled for the first few pictures.

For generating the random points in the pixel and the light source I implemented stratified sampling. I calculated the total area to sample and the parameters defining that part of the view in order to evenly divide it into an n x n grid of patches, where n = ceiling(sqrt(number of rays)). I created a vector of grid points and randomly chose from them to send out a ray in a random spot inside that area. I did this for as many times as there were rays to send out.



RADIOSITY:
My form factor implementation seems to be correct, but falls short in the final values. The form factor was calculated using an approximation of surface to surface via one point. This was fixed by scaling the form factors for a face with all other faces to always sum to 1.

The cornell box scene appears to properly render. I did not notice any differences between what I got and what the homework page showed. The radiosity solution was an iterative process.

For visibility/occlusion, I was able to get something working, though it is not complete. The lighting is very dim when it is rendered, and the back area is completely dark. Subdividing further seems to help improve the rendering.

The order notation of visibility/occlusion testing is O(n^3), as n^2 comparisons need to be made between the patches with n raycasts for shadows.

PHOTON MAPPING:

-Did not implement


OTHER NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.

Stratified sampling:
My sampling implementation took only a short while longer to render the reflective sphere and textured plane scene. Though there was no noticeable quality increase with and without on only soft shadows with 4 samples.

With 9 shadow samples and 9 alias samples, stratified sampling took 10 minutes, while without it took 6. With the sampling on the antialiasing the quality was significantly better.

My code is set to only use stratified sampling. To run without it requires commenting out sections.

Glossy surfaces: 
Input command: 
-input ../src/textured_plane_reflective_sphere.obj -num_bounces 1 -num_shadow_samples 9 -num_antialias_samples 9 -num_glossy_samples 9

I implemented glossy surfaces by casting additional rays when casting a reflection ray. These additional rays are cast in a similar direction, but are moved slightly to sample the area around the reflection ray. I ran my renderer with antialiasing and soft shadows to produce the best result possible. With the specified command line my renderer took just under 20 minutes to finish.