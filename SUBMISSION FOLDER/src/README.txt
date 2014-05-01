HOMEWORK 3: RAY TRACING, RADIOSITY, & PHOTON MAPPING

NAME:  < Gerrett Diamond >

Version: GLUT
OS: Windows

TOTAL TIME SPENT:  < 30 >
Please estimate the number of hours you spent on this assignment.


COLLABORATORS: 
You must do this assignment on your own, as described in the 
Academic Integrity Policy.  If you did discuss the problem or errors
messages, etc. with anyone, please list their names here.



RAYTRACING:
< insert comments on the implementation, known bugs, extra credit >
***Random point generation***
For the light source I just used the random point function that was already built within the Face class to draw the shadow segements. The pixel was done by generating a random section of the pixel using the Mersenne Twister generator.


RADIOSITY:
< insert comments on the implementation, known bugs, extra credit >
***Form factor computation***
	So the form factor calculation ends up with an order notation of O(n^3*f*s) where n is the number of faces, f the number of form factor samples, and s the number of shadow samples. n is cubed due to the comparison of each face to each face and then the casyray also includes a third n. 
	Increasing the number of faces, although it slows down drastically is important to give smoother boundaries between the color bleeding. For just the cornell_box.obj test case, if keep at the original faces, the walls become a solid purple color but after subdividing a few times it can be seen that the color bleeding fades across the wall in a smoother fashion. 


PHOTON MAPPING:
< insert comments on the implementation, known bugs, extra credit >
I only gave myself a 1 on the ray tracing with photons because I didnt think that I was getting the correct behavior with the extra lighting.
Also the program occasionally crashes for larger photon counts but then works on other times.


OTHER NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.

I added a test case mytest.obj for ray tracing. It consists of 3 colored spheres (red, green, and blue) and a black sphere sitting on top of a plane of light sources each with different coloring which was made distiguishable as red, blue, green, and white. It takes a high number of samples to get good images and a some time so some images were taken (mytestside.png, mytesttop.png).
./render -input ../src/mytest.obj -num_bounces 5 -num_shadow_samples 10 -num_antialias_samples 20
