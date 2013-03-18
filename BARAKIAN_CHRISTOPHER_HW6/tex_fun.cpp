/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"math.h"

GzColor	*image;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;
  float s,t;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
  if( u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0){
	  color[RED] = 0.0;
	  color[GREEN] = 0.0;
	  color[BLUE] = 0.0;
	  return 0;
  }
  else{
	  /* determine texture cell corner values and perform bilinear interpolation */
	  int minX = u*(xs-1);
	  int minY = v*(ys-1);
	  int maxX = minX + 1;
	  int maxY = minY + 1;
	  s = u*(xs-1) - minX;
	  t = v*(ys-1) - minY;
	  int a = (xs*minY)+minX;
	  int b = (xs*minY)+maxX;
	  int c = (xs*maxY)+maxX;
	  int d = (xs*maxY)+minX;

	  /* set color to interpolated GzColor value and return */
	  color[RED] = (s*t*image[c][RED]) + ((1-s)*t*image[d][RED]) + (s*(1-t)*image[b][RED]) + ((1-s)*(1-t)*image[a][RED]);
	  color[GREEN] = (s*t*image[c][GREEN]) + ((1-s)*t*image[d][GREEN]) + (s*(1-t)*image[b][GREEN]) + ((1-s)*(1-t)*image[a][GREEN]);
	  color[BLUE] = (s*t*image[c][BLUE]) + ((1-s)*t*image[d][BLUE]) + (s*(1-t)*image[b][BLUE]) + ((1-s)*(1-t)*image[a][BLUE]);

	  return 0;
  }


}


/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	/*if(u < .5 && v < .5){
		color[0] = 0;
		color[1] = 0;
		color[2] = 0;
	}
	else if(u > .5 && v < .5){
		color[0] = 1;
		color[1] = 1;
		color[2] = 1;
	}
	else if(u < .5 && v > .5){
		color[0] = 1;
		color[1] = 1;
		color[2] = 1;
	}
	else{
		color[0] = 0;
		color[1] = 0;
		color[2] = 0;
	}*/

	int n = u*1;
    n = (n<<13) ^ n;
    float noise1 = ( 1.0 - ( (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 1073741824.0);
	n = (u*256)+1;
    n = (n<<13) ^ n;
    float noise2 = ( 1.0 - ( (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 1073741824.0);

	float weight = (u*1) - floor(u*1);
	float noise = noise1*(1-weight)+(noise2*weight);

	n = v*1;
    n = (n<<13) ^ n;
    float noise3 = ( 1.0 - ( (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 1073741824.0);
	n = (v*256)+1;
    n = (n<<13) ^ n;
    float noise4 = ( 1.0 - ( (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 1073741824.0);

	weight = (v*1) - floor(v*1);
	float noisev = noise3*(1-weight)+(noise4*weight);

	/*color[0] = sin(3.14*(noise)/noisev)+cos(3.14*(noise)/noisev);
	color[1] = sin(3.14*noise)+cos(3.14*noisev);
	color[2] = cos(3.14*(noise)/noisev);*/

	color[0] = (fmod(floor(u*16)+floor(v*16),2) < 1) ? sin(3.14*(noise)/noisev)+cos(3.14*(noise)/noisev):cos(3.14*v)/4;
	color[1] = (fmod(floor(u*16)+floor(v*16),2) < 1) ? sin(3.14*noise)+cos(3.14*noisev):cos(3.14*u*2);
	color[2] = (fmod(floor(u*16)+floor(v*16),2) < 1) ? cos(3.14*(noise)/noisev):sin(3.14*v*2);

	return 0;
}

