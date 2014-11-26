/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	"math.h"
GzColor	*image=NULL;
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

  if (reset) {          /* open and load texture file */
    fd = fopen ("WindowsImage", "rb");
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
  if (u > 1)
  {
	  u = 1;
  }
  if (u < 0)
  {
	  u = 0;
  }
  if (v > 1)
  {
	  v = 1;
  }
  if (v < 0)
  {
	  v = 0;
  }
/* determine texture cell corner values and perform bilinear interpolation */
  float px, py, s, t;
  px = u * (xs - 1);
  py = v * (ys - 1);
  s = px - floor(px);
  t = py - floor(py);

/* set color to interpolated GzColor value and return */
  GzColor colorA;
  GzColor colorB;
  GzColor colorC;
  GzColor colorD;

  for (i = 0; i < 3; i++)
  {
	  colorA[i] = image[xs * (int)floor(py) + (int)floor(px)][i];
  }

  for (i = 0; i < 3; i++)
  {
	  colorB[i] = image[xs * (int)floor(py) + (int)ceil(px)][i];
  }

  for (i = 0; i < 3; i++)
  {
	  colorC[i] = image[xs * (int)ceil(py) + (int)ceil(px)][i];
  }

  for (i = 0; i < 3; i++)
  {
	  colorD[i] = image[xs * (int)ceil(py) + (int)floor(px)][i];
  }

  for (i = 0; i < 3; i++)
  {
	  color[i] = s * t * colorC[i] + (1 - s) * t * colorD[i] + s * (1 - t) * colorB[i] + (1 - s) * (1 - t) * colorA[i]; 
  }
  return GZ_SUCCESS;
}
struct complexNumber
{
public:
float r;
float i;
};

/* Procedural texture function */
//Fractal texture: Julia Set
int ptex_fun(float u, float v, GzColor color)
{
	u = u*255;
	v = v*255;	

	int Box_X = 16;
	if((int)((int)u/Box_X + (int)v/Box_X)%2 == 0)
	{
		color[0] = 0.3;
		color[1] = 0.4;
		color[2] = 0.5;
	}
	else
	{
		color[0] = 0.5;
		color[1] = 0.6;
		color[2] = 0.7;
	}
    return GZ_SUCCESS;
}
/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);

	return GZ_SUCCESS;
}
