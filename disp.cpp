/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"

int CheckForBounds(int _Value, int _LowerBounds, int _UpperBounds)
{
	if(_Value < _LowerBounds)
		_Value = _LowerBounds;
	if(_Value > _UpperBounds)
		_Value = _UpperBounds;

	return _Value;
}

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
	//ALLOCATING MEMORY FOR FRAME BUFFER
	*framebuffer = (char*) malloc(width * height * sizeof(GzPixel));

	return GZ_SUCCESS;
}

int   GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
	//CREATING A NEW DISPLAY AND ALLOCATIN MEMORY FOR FRAME BUFFER
	*display = (GzDisplay*) malloc(sizeof(GzDisplay));
	(*display)->fbuf = (GzPixel*) malloc(xRes * xRes * sizeof(GzPixel));
	(*display)->xres = xRes;
	(*display)->yres = yRes;

	return GZ_SUCCESS;
}
int GzFreeDisplay(GzDisplay	*display)
{
	free(display->fbuf);
	free(display);
	return GZ_SUCCESS;
}
int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
	xRes = (int*)display->xres;
	yRes = (int*)display->yres;

	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
	//INITIALIZING THE DISPLAY TO BLACK PIXELS
	GzPixel *_TempBuffer = display->fbuf;
	for(int i=0 ; i<display->xres ; i++)
	{
		for(int j=0 ; j<display->xres ; j++)
		{
			_TempBuffer = display->fbuf + ARRAY(j,i);
			_TempBuffer->alpha = 1;
			_TempBuffer->red = 0;
			_TempBuffer->blue = 000;
			_TempBuffer->green = 1000;
			_TempBuffer->z = INT_MAX;
		}
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	//CHECKING FOR BOUNDS
	r = CheckForBounds(r, 0, 4095);
	g = CheckForBounds(g, 0, 4095);
	b = CheckForBounds(b, 0, 4095);

	i = CheckForBounds(i, 0, display->xres-1);
	j = CheckForBounds(j, 0, display->yres-1);

	//ASSIGNING THE RGBAZ VALUES
	GzPixel *tempBuffer = display->fbuf;
	tempBuffer = display->fbuf + ARRAY(i,j);
	tempBuffer->alpha = 1;
	tempBuffer->red = r;
	tempBuffer->blue = b;
	tempBuffer->green = g;
	tempBuffer->z = z;

	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	if(i > display->xres || j > display->yres || i < 0 || j < 0)
		return GZ_FAILURE;

	GzPixel *tempBuffer = display->fbuf + ARRAY(i,j);
	*r = tempBuffer->red;
	*g = tempBuffer->green;
	*b = tempBuffer->blue;
	*a = tempBuffer->alpha;
	*z = tempBuffer->z;

	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
	//WRITING THE VALUES TO FILE
	fprintf(outfile,"P6 %d %d 255\r",display->xres, display->yres);

	GzPixel* _StartPointer = display->fbuf;
	for(int _RowNo=0 ; _RowNo<display->xres; _RowNo++)
	{
		for(int _ColumnNo=0; _ColumnNo<display->yres ; _ColumnNo++)
		{
			_StartPointer = display->fbuf + ARRAY(_ColumnNo,_RowNo);
			fprintf(outfile, "%c%c%c", _StartPointer->red>>4, _StartPointer->green>>4, _StartPointer->blue>>4);
		}
	} 
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{
	//WRITING VALUES TO FRAME BUFFER
	GzPixel *_TempBuffer = display->fbuf;
	for(int _RowNo =0 ; _RowNo<display->xres ; _RowNo++)
	{
		for(int _ColumnNo=0 ; _ColumnNo<display->yres ; _ColumnNo++)
		{
			_TempBuffer = display->fbuf + ARRAY(_ColumnNo,_RowNo);
			*(framebuffer++) = _TempBuffer->blue>>4;
			*(framebuffer++) = _TempBuffer->green>>4;
			*(framebuffer++) = _TempBuffer->red>>4;
		}
	}
	return GZ_SUCCESS;
}