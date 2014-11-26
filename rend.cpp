/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"


typedef struct point
{
	float x,y,z;
	void SetValues(float _x,float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

}Point;

typedef struct Normals
{
	float x,y,z;
	void SetValues(float _x,float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
}Normals;

struct Vector3
{
	float m_x;
	float m_y;
	float m_z;
	float nx, ny, nz;
	Vector3() : m_x(0), m_y(0), m_z(0){}

	Vector3(float _x, float _y, float _z)
	{
		m_x = _x;
		m_y = _y;
		m_z = _z;
	}
	Vector3 Subtract(const Vector3 &v)
	{
		Vector3 _Vec = Vector3(m_x -v.m_x, m_y- v.m_y, m_z -v.m_x);
		return 	_Vec;

	}
	Vector3 Add(const Vector3 &v)
	{
		Vector3 _Vec = Vector3(m_x +v.m_x, m_y+ v.m_y, m_z +v.m_x);
		return 	_Vec;

	}
	Vector3 MultiplyByScalor(float value)
	{
		Vector3 _Vec = Vector3(m_x * value, m_y * value, m_z * value);
		return 	_Vec;

	}
	Vector3 Cross(const Vector3 &v)
	{
		Vector3 _Vec = Vector3(m_y * v.m_z - m_z * v.m_y, m_z * v.m_x - m_x * v.m_z, m_x * v.m_y - m_y * v.m_x);
		return 	_Vec;
	}

	inline Vector3 operator / ( float f) {return Vector3(m_x / f, m_y / f, m_z / f);}
	void operator /= (const float &f) {m_x /= f; m_y /= f; m_z /= f;}

	inline Vector3 operator * ( float f) {return Vector3(m_x * f, m_y * f, m_z * f);}
	void operator *=(const float &f) {m_x *= f; m_y *= f; m_z *= f;}

	float GetLength()
	{
		return sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
	}
	float Dot(const Vector3 &v) 
	{ 
		return m_x * v.m_x + m_y * v.m_y + m_z * v.m_z;
	}
	void Normalise()
	{
		float length = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
		m_x = m_x / length;
		m_y = m_y / length;
		m_z = m_z / length;
	}
};
struct Plane
{
	float m_A,m_B,m_C,m_D;
	float m_Ax,m_Ay,m_Az;
	float m_Bx,m_By,m_Bz;
	float m_AC, m_BC;

	Vector3 m_Normal;
	Plane() : m_A(0), m_B(0), m_C(0), m_D(0){}
	void CalculateValues(Point _VertexPosition[3])
	{
		m_Ax = _VertexPosition[1].x - _VertexPosition[0].x;
		m_Ay = _VertexPosition[1].y - _VertexPosition[0].y;
		m_Az = _VertexPosition[1].z - _VertexPosition[0].z;
		m_Bx = _VertexPosition[2].x - _VertexPosition[1].x;
		m_By = _VertexPosition[2].y - _VertexPosition[1].y;
		m_Bz = _VertexPosition[2].z - _VertexPosition[1].z;

		m_A = m_Ay * m_Bz - m_Az * m_By ;
		m_B = m_Az * m_Bx - m_Ax * m_Bz ;
		m_C = m_Ax * m_By - m_Ay * m_Bx ;
		m_D = - m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_VertexPosition[0].z) ;
	}
	void CalculateColor(Point _VertexPosition[3], GzColor _ColorVertex[3], int _Value)
	{
		m_Ax = _VertexPosition[1].x - _VertexPosition[0].x;
		m_Ay = _VertexPosition[1].y - _VertexPosition[0].y;
		m_AC = _ColorVertex[1][_Value] - _ColorVertex[0][_Value];
		m_Bx = _VertexPosition[2].x - _VertexPosition[1].x;
		m_By = _VertexPosition[2].y - _VertexPosition[1].y;
		m_BC = _ColorVertex[2][_Value] - _ColorVertex[1][_Value];
			
		m_A = m_Ay * m_BC - m_AC * m_By ;
		m_B = m_AC * m_Bx - m_Ax * m_BC ;
		m_C = m_Ax * m_By - m_Ay * m_Bx ;
		m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_ColorVertex[0][_Value]) ;
	}
	void CalculateNormals(Point _VertexPosition[3], Normals _ColorVertex[3], int _Value)
	{
		m_Ax = _VertexPosition[1].x - _VertexPosition[0].x;
		m_Ay = _VertexPosition[1].y - _VertexPosition[0].y;
	
		m_Bx = _VertexPosition[2].x - _VertexPosition[1].x;
		m_By = _VertexPosition[2].y - _VertexPosition[1].y;

		if(_Value == 0)
		{
			m_AC = _ColorVertex[1].x - _ColorVertex[0].x;
			m_BC = _ColorVertex[2].x - _ColorVertex[1].x;
		}
		else if(_Value == 1)
		{
			m_AC = _ColorVertex[1].y - _ColorVertex[0].y;
			m_BC = _ColorVertex[2].y - _ColorVertex[1].y;
		}
		else if(_Value == 2)
		{
			m_AC = _ColorVertex[1].z - _ColorVertex[0].z;
			m_BC = _ColorVertex[2].z - _ColorVertex[1].z;
		}

			
		m_A = m_Ay * m_BC - m_AC * m_By ;
		m_B = m_AC * m_Bx - m_Ax * m_BC ;
		m_C = m_Ax * m_By - m_Ay * m_Bx ;
		if(_Value == 0)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_ColorVertex[0].x) ;
		else if(_Value == 1)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_ColorVertex[0].y) ;
		else if(_Value == 2)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_ColorVertex[0].z) ;
	}
	void CalculateTextures(Point _VertexPosition[3], Point _TextureVertex[3], int _Value)
	{
		m_Ax = _VertexPosition[1].x - _VertexPosition[0].x;
		m_Ay = _VertexPosition[1].y - _VertexPosition[0].y;
	
		m_Bx = _VertexPosition[2].x - _VertexPosition[1].x;
		m_By = _VertexPosition[2].y - _VertexPosition[1].y;

		if(_Value == 0)
		{
			m_AC = _TextureVertex[1].x - _TextureVertex[0].x;
			m_BC = _TextureVertex[2].x - _TextureVertex[1].x;
		}
		else if(_Value == 1)
		{
			m_AC = _TextureVertex[1].y - _TextureVertex[0].y;
			m_BC = _TextureVertex[2].y - _TextureVertex[1].y;
		}
		else if(_Value == 2)
		{
			m_AC = _TextureVertex[1].z - _TextureVertex[0].z;
			m_BC = _TextureVertex[2].z - _TextureVertex[1].z;
		}

			
		m_A = m_Ay * m_BC - m_AC * m_By ;
		m_B = m_AC * m_Bx - m_Ax * m_BC ;
		m_C = m_Ax * m_By - m_Ay * m_Bx ;
		if(_Value == 0)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_TextureVertex[0].x) ;
		else if(_Value == 1)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_TextureVertex[0].y) ;
		else if(_Value == 2)
			m_D = -m_A * (_VertexPosition[0].x) - m_B * (_VertexPosition[0].y) - m_C * (_TextureVertex[0].z) ;
	}
	void CalculateAndNormalize()
	{
		float mag = sqrt(m_A * m_A + m_B * m_B + m_C * m_C);
		m_Normal = Vector3(m_A/mag, m_B/mag, m_C/mag);
	
		m_Normal.Normalise();
	}
};

void RefractedRayCalculation(GzRender *render, Normals _NormalTemp, Vector3 *refractedRay)
{
	Vector3 _NormalVector = Vector3(_NormalTemp.x, _NormalTemp.y, _NormalTemp.z);
	Vector3 _Eye = Vector3(0, 0, -1);
	float n = 1.5;

	// k = 1.0 - eta * eta * (1.0 - dot(N, I) * dot(N, I));
	float VdotN =  _Eye.Dot(_NormalVector);
	float k = 1.0 - n * n * (1.0 - VdotN * VdotN);

	// R = eta * I - (eta * dot(N, I) + sqrt(k)) * N;
	Vector3  RefractedRay;
	RefractedRay = _Eye.MultiplyByScalor(n).Subtract(_NormalVector.MultiplyByScalor(n * VdotN + sqrt(k)));

	refractedRay->m_x = RefractedRay.m_x;
	refractedRay->m_y = RefractedRay.m_y;
	refractedRay->m_z = RefractedRay.m_z;
}

// old one
	/*
	Vector3 temp = Vector3(_NormalVector.m_x*VdotN, _NormalVector.m_y*VdotN, _NormalVector.m_z*VdotN); //N* NdotE
	temp = _Eye.Subtract(temp);

	float val = sqrt(k);
	Vector3 temp1 = _NormalVector;
	temp1.m_x *= val;
	temp1.m_y *= val;
	temp1.m_z *= val;

	Vector3  RefractedRay = (temp*n).Add(temp1);
	*/
	//sprintf(debug_buf,"\n %f %f %f", RefractedRay.m_x, RefractedRay.m_y,RefractedRay.m_z);
	//OutputDebugString(debug_buf);

int normalMatrixlevel = 0;
void GzMatrixMutiply(GzMatrix &result,GzMatrix first , GzMatrix second)
{
	float _Sum = 0;	 
	for ( int i = 0 ; i < 4 ; i++ )
	{
		for ( int j = 0 ; j < 4 ; j++ )
		{
			for ( int k = 0 ; k < 4 ; k++ )
		    {
				_Sum = _Sum + first[i][k]*second[k][j];
		    }
		    result[i][j] = _Sum;
		    _Sum = 0;
		}
	}
}


/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}
int VectorNormalise(GzCoord vector)
{
	if (vector == NULL)
	{
		return GZ_FAILURE;
	}

	float length = sqrt((vector[X]*vector[X]) + (vector[Y]*vector[Y]) + (vector[Z]*vector[Z]));

	vector[X] = vector[X] / length;
	vector[Y] = vector[Y] / length;
	vector[Z] = vector[Z] / length;

	return GZ_SUCCESS;
}
void InitializeMatrix(GzMatrix *_Matrix)
{
	for(int i =0; i< 4; i++)
	{
		for(int j=0; j<4; j++)
		{
			(*_Matrix)[i][j] =0.0;
		}
	}
}
void GetIdentitiyMatrix(GzMatrix* _Identity)
{
	(*_Identity)[0][0] = 1;	
	(*_Identity)[0][1] = 0;	
	(*_Identity)[0][2] = 0;	    
	(*_Identity)[0][3] = 0;

	(*_Identity)[1][0] = 0;	
	(*_Identity)[1][1] = 1;	
	(*_Identity)[1][2] = 0;	   	
	(*_Identity)[1][3] = 0;

	(*_Identity)[2][0] = 0;	
	(*_Identity)[2][1] = 0;	
	(*_Identity)[2][2] = 1;		
	(*_Identity)[2][3] = 0;
	
	(*_Identity)[3][0] = 0;	
	(*_Identity)[3][1] = 0;	
	(*_Identity)[3][2] = 0;		
	(*_Identity)[3][3] = 1;
}
int GzRotXMat(float _Degree, GzMatrix _Matrix)
{
	InitializeMatrix((GzMatrix*)_Matrix);
	
	double rad = _Degree * 3.142 / 180.0;
	_Matrix[0][0] = 1.0;
	_Matrix[1][1] = cos(rad);
	_Matrix[1][2] = -sin(rad);
	_Matrix[2][1] = sin(rad);
	_Matrix[2][2] = cos(rad);
	_Matrix[3][3] = 1.0;
	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix _Matrix)
{
	InitializeMatrix((GzMatrix*)_Matrix);

	double rad = degree * 3.142 / 180.0;
	_Matrix[0][0] = cos(rad);
	_Matrix[0][2] = sin(rad);
	_Matrix[1][1] = 1.0;
	_Matrix[2][0] = -sin(rad);
	_Matrix[2][2] = cos(rad);
	_Matrix[3][3] = 1.0;

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix _Matrix)
{
	InitializeMatrix((GzMatrix*)_Matrix);

	double rad = degree * 3.142 / 180.0;
	_Matrix[0][0] = cos(rad);
	_Matrix[0][1] = -sin(rad);
	_Matrix[1][0] = sin(rad);
	_Matrix[1][1] = cos(rad);
	_Matrix[2][2] = 1.0;
	_Matrix[3][3] = 1.0;

	return GZ_SUCCESS;
}



int GzTrxMat(GzCoord translate, GzMatrix _Matrix)
{
	InitializeMatrix((GzMatrix*)_Matrix);

	_Matrix[0][0] =	1;			    _Matrix[0][3] = translate[0];
	_Matrix[1][1] =	1;			    _Matrix[1][3] = translate[1];
	_Matrix[2][2] = 1;				_Matrix[2][3] = translate[2];
	_Matrix[3][3] = 1;

	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix _Matrix)
{
	InitializeMatrix((GzMatrix*)_Matrix);

	_Matrix[0][0] =	scale[0];		
	_Matrix[1][1] =	scale[1];		
	_Matrix[2][2] = scale[2];		
	_Matrix[3][3] = 1;

	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render,  GzDisplay	*display)
{
	//**************CREATING INITIAL VALUES
	GzRender *rend;
	rend = new GzRender;
	rend->display = display;
	if(!rend)
		return GZ_FAILURE;
	
	
	rend->numlights = 0;
	rend->matlevel= -1;
	
	rend->camera.FOV = DEFAULT_FOV;
	rend->camera.lookat[0] = 0.0;
	rend->camera.lookat[1] = 0.0;
	rend->camera.lookat[2] = 0.0;
	rend->camera.position[0] = DEFAULT_IM_X;
	rend->camera.position[1] = DEFAULT_IM_Y;
	rend->camera.position[2] = DEFAULT_IM_Z;
	rend->camera.worldup[0] = 0;
	rend->camera.worldup[1] = 1;
	rend->camera.worldup[2] = 0;

		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		{
			rend->Xsp[i][j] = 0;
			rend->camera.Xiw[i][j] = 0;
			rend->camera.Xpi[i][j] = 0;
		}
	//GzCamera *camera;
	//camera = new GzCamera;
	//GzPutCamera(rend,camera);

	*render = rend;
	return GZ_SUCCESS;
}
void ClampColors(GzColor *_Color)
{
}
int GzFreeRender(GzRender *render)
{
	free(render);

	return GZ_SUCCESS;
}

float CalculateMag(GzCoord vec)
{
	return sqrtf((vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

void shade(GzCoord norm, GzCoord color)
{
  GzCoord	light;
  float		coef;

  light[0] = 0.707f;
  light[1] = 0.5f;
  light[2] = 0.5f;

  coef = light[0]*norm[0] + light[1]*norm[1] + light[2]*norm[2];
  if (coef < 0) 	coef *= -1;

  if (coef > 1.0)	coef = 1.0;
  color[0] = coef*0.95f;
  color[1] = coef*0.65f;
  color[2] = coef*0.88f;
}

int GzBeginRender(GzRender *render)
{
	GzInitDisplay(render->display);
	float _Radians = render->camera.FOV * (3.1414 / 180);
	float _InverseDistance = tan (_Radians/2);

	//CALCULATING XSP 
	render->Xsp[0][0] = render->display->xres/2;
	render->Xsp[0][3] = render->display->xres/2;
	render->Xsp[1][1] = -render->display->yres/2;
	render->Xsp[1][3] = render->display->yres/2;
	render->Xsp[2][2] = INT_MAX * tan((render->camera.FOV / 2.0) * (3.142 / 180.0));
	render->Xsp[3][3] = 1.0;


	//CALCULATING THE INITIAL VALUES OF CAMERA
	render->camera.Xpi[0][0] = 1;		
	render->camera.Xpi[0][1] = 0;		
	render->camera.Xpi[0][2] = 0;				
	render->camera.Xpi[0][3] = 0;
	render->camera.Xpi[1][0] = 0;		
	render->camera.Xpi[1][1] = 1;		
	render->camera.Xpi[1][2] = 0;	   			
	render->camera.Xpi[1][3] = 0;
	render->camera.Xpi[2][0] = 0;		
	render->camera.Xpi[2][1] = 0;		
	render->camera.Xpi[2][2] = 1;	
	render->camera.Xpi[2][3] = 0;
	render->camera.Xpi[3][0] = 0;		
	render->camera.Xpi[3][1] = 0;		
	render->camera.Xpi[3][2] = _InverseDistance;	
	render->camera.Xpi[3][3] = 1;

	GzCoord _AxisX;
	GzCoord	_AxisY;
	GzCoord _AxisZ;

	_AxisZ[X] = render->camera.lookat[X] - render->camera.position[X];
	_AxisZ[Y] = render->camera.lookat[Y] - render->camera.position[Y];
	_AxisZ[Z] = render->camera.lookat[Z] - render->camera.position[Z];

	float mod = CalculateMag(_AxisZ);

	if(CalculateMag(_AxisZ) != 0)
	{
		_AxisZ[X] = _AxisZ[X] / mod;
		_AxisZ[Y] = _AxisZ[Y] / mod;
		_AxisZ[Z] = _AxisZ[Z] / mod;
	}
	
	GzCoord _NewUp,_Vector1, _Vector2;
	float dot_product;

	dot_product = render->camera.worldup[X] * _AxisZ[X] +  render->camera.worldup[Y] * _AxisZ[Y] +  render->camera.worldup[Z] * _AxisZ[Z];

	_Vector2[X] = _AxisZ[X] * dot_product;
	_Vector2[Y] = _AxisZ[Y] * dot_product;
	_Vector2[Z] = _AxisZ[Z] * dot_product;

	_NewUp[X] = render->camera.worldup[X] - _Vector2[X];
	_NewUp[Y] = render->camera.worldup[Y] - _Vector2[Y];
	_NewUp[Z] = render->camera.worldup[Z] - _Vector2[Z];

	mod = CalculateMag(_NewUp);

    _AxisY[X] = _NewUp[X] / mod;
	_AxisY[Y] = _NewUp[Y] / mod;
	_AxisY[Z] = _NewUp[Z] / mod;

	_AxisX[X] = _AxisY[1] * _AxisZ[2] - _AxisY[2] * _AxisZ[1];
	_AxisX[Y] = _AxisY[2] * _AxisZ[0] - _AxisY[0] * _AxisZ[2];
	_AxisX[Z] = _AxisY[0] * _AxisZ[1] - _AxisY[1] * _AxisZ[0];

	VectorNormalise(_AxisX);
	VectorNormalise(_AxisY);
	VectorNormalise(_AxisZ);

	//CALCULATING THE DOT PRODUCTS WITHT HE CAMERA
	float _DotXC = _AxisX[X] * render->camera.position[X] + _AxisX[Y] * render->camera.position[Y] + _AxisX[Z] * render->camera.position[Z];
	float _DotYC = _AxisY[X] * render->camera.position[X] + _AxisY[Y] * render->camera.position[Y] + _AxisY[Z] * render->camera.position[Z];
	float _DotZC = _AxisZ[X] * render->camera.position[X] + _AxisZ[Y] * render->camera.position[Y] + _AxisZ[Z] * render->camera.position[Z];

	render->camera.Xiw[0][0] =	_AxisX[0];	
	render->camera.Xiw[0][1] =	_AxisX[1];		
	render->camera.Xiw[0][2] =	_AxisX[2];	    
	render->camera.Xiw[0][3] = -_DotXC;

	render->camera.Xiw[1][0] =	_AxisY[0];	
	render->camera.Xiw[1][1] =	_AxisY[1];		
	render->camera.Xiw[1][2] =  _AxisY[2];	   	
	render->camera.Xiw[1][3] = -_DotYC;

	render->camera.Xiw[2][0] =	_AxisZ[0];	
	render->camera.Xiw[2][1] =	_AxisZ[1];		
	render->camera.Xiw[2][2] =  _AxisZ[2];		
	render->camera.Xiw[2][3] = -_DotZC;

	render->camera.Xiw[3][0] = 0;	
	render->camera.Xiw[3][1] = 0;		
	render->camera.Xiw[3][2] = 0;		
	render->camera.Xiw[3][3] = 1;

	GzMatrix _Xnormiw;

	_Xnormiw[0][0] = _AxisX[X];	
	_Xnormiw[0][1] = _AxisX[Y];		
	_Xnormiw[0][2] = _AxisX[Z];	    
	_Xnormiw[0][3] = 0;

	_Xnormiw[1][0] = _AxisY[X];	
	_Xnormiw[1][1] = _AxisY[Y];		
	_Xnormiw[1][2] = _AxisY[Z];	   	
	_Xnormiw[1][3] = 0;

	_Xnormiw[2][0] = _AxisZ[X];	
	_Xnormiw[2][1] = _AxisZ[Y];		
	_Xnormiw[2][2] = _AxisZ[Z];		
	_Xnormiw[2][3] = 0;

	_Xnormiw[3][0] = 0;	
	_Xnormiw[3][1] = 0;		
	_Xnormiw[3][2] = 0;		
	_Xnormiw[3][3] = 1;


	//STACK WITH IDENTITY  MATRIX
	GzMatrix _Identity;
	GetIdentitiyMatrix((GzMatrix*)_Identity);
	render->matlevel = -1;
	
	GzMatrix _TempMatrix;

	render->matlevel++;
	memcpy(render->Ximage[render->matlevel],_Identity,sizeof(GzMatrix));
	
	normalMatrixlevel = 0;
	memcpy(render->Xnorm[normalMatrixlevel],_Identity,sizeof(GzMatrix));

	//Push Xsp
	render->matlevel++;
	GzMatrixMutiply(_TempMatrix,render->Ximage[render->matlevel-1],render->Xsp);
	memcpy(render->Ximage[render->matlevel],_TempMatrix,sizeof(GzMatrix));

	//Push Xpi
	render->matlevel++;
	GzMatrixMutiply(_TempMatrix,render->Ximage[render->matlevel-1],render->camera.Xpi);
	memcpy(render->Ximage[render->matlevel],_TempMatrix,sizeof(GzMatrix));

	//Push Xiw
	render->matlevel++;
	GzMatrixMutiply(_TempMatrix,render->Ximage[render->matlevel-1],render->camera.Xiw);
	memcpy(render->Ximage[render->matlevel],_TempMatrix,sizeof(GzMatrix));


	GzMatrixMutiply(_TempMatrix,render->Xnorm[normalMatrixlevel],_Xnormiw);
	memcpy(render->Xnorm[++normalMatrixlevel],_TempMatrix,sizeof(GzMatrix));

	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition

*/
	if (render == NULL) {
		return GZ_FAILURE;
	}
	if (camera == NULL) {
		return GZ_FAILURE;
	}

	render->camera.FOV = camera->FOV;
	render->camera.lookat[X] = camera->lookat[X];
	render->camera.lookat[Y] = camera->lookat[Y];
	render->camera.lookat[Z] = camera->lookat[Z];
	render->camera.position[X] = camera->position[X];
	render->camera.position[Y] = camera->position[Y];
	render->camera.position[Z] = camera->position[Z];
	render->camera.worldup[X] = camera->worldup[X]; 
	render->camera.worldup[Y] = camera->worldup[Y]; 
	render->camera.worldup[Z] = camera->worldup[Z]; 

	return GZ_SUCCESS;	
}

int NormalizeMatrix(GzMatrix matrix)
{
	float a,b,c;
	
	for(int i=0;i<3;i++)
	{

		matrix[i][3] = 0;

		a = matrix[i][0];
		b = matrix[i][1];
		c = matrix[i][2];

		float normfactor = 1 / (sqrtf( a * a + b * b + c * c) );

		matrix[i][0] = normfactor * matrix[i][0];
		matrix[i][1] = normfactor * matrix[i][1];
		matrix[i][2] = normfactor * matrix[i][2];

	}

	return GZ_SUCCESS;
}

int GzPushNormalMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/

	if (render == NULL) {
		return GZ_FAILURE;
	}
	if (matrix == NULL) {
		return GZ_FAILURE;
	}

	GzMatrix temp;

	NormalizeMatrix(matrix);

    GzMatrixMutiply(temp,render->Xnorm[normalMatrixlevel],matrix);
	memcpy(render->Xnorm[++normalMatrixlevel],temp,sizeof(GzMatrix));

	NormalizeMatrix(matrix);

	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/

	if (render == NULL) {
		return GZ_FAILURE;
	}
	if (matrix == NULL) {
		return GZ_FAILURE;
	}

	GzMatrix temp;
	render->matlevel++;

    GzMatrixMutiply(temp,render->Ximage[render->matlevel-1],matrix);
	memcpy(render->Ximage[render->matlevel],temp,sizeof(GzMatrix));

	GzPushNormalMatrix(render,matrix);

	return GZ_SUCCESS;
}
int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/

	if (render == NULL) {
		return GZ_FAILURE;
	}
	
	if (render->Ximage[0] == NULL) {
		return GZ_FAILURE;
	}
	if (render->matlevel < 0) {
		return GZ_FAILURE;
	}

	render->matlevel -=1;

	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	if (render == NULL)  
	{
		return GZ_FAILURE;
	}
	if (nameList == NULL) 
	{
		return GZ_FAILURE;
	}
	if (valueList == NULL) 
	{
		return GZ_FAILURE;
	}

	for (int i = 0; i < numAttributes; i++)
	{
		
		if (nameList[i] == GZ_DIRECTIONAL_LIGHT)
		{
			render->lights[render->numlights] = *(GzLight*)(valueList[i]);
			render->numlights++ ;
		}
		
		else if (nameList[i] == GZ_AMBIENT_LIGHT)
		{
			render->ambientlight = *(GzLight*)(valueList[i]);
		}
		
		
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT)
		{
			GzColor* _DiffuseColor = (GzColor*)(valueList[i]);
			render->Kd[0] = _DiffuseColor[0][0];
			render->Kd[1] = _DiffuseColor[0][1];
			render->Kd[2] = _DiffuseColor[0][2];
		}

		else if (nameList[i] == GZ_INTERPOLATE)
		{
			render->interp_mode = *(int*)valueList[i];
		}
		
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT)
		{
			GzColor* _AmbientColor = (GzColor*)(valueList[i]);
			render->Ka[0] = _AmbientColor[0][0];
			render->Ka[1] = _AmbientColor[0][1];
			render->Ka[2] = _AmbientColor[0][2];
		}

		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT)
		{
			GzColor* specColor = (GzColor*)(valueList[i]);
			render->Ks[0] = specColor[0][0];
			render->Ks[1] = specColor[0][1];
			render->Ks[2] = specColor[0][2];
		}

		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT)
		{
			render->spec = *(float*)valueList[i];
		}

		
		//render->Kd

		else if (nameList[i] == GZ_RGB_COLOR)
		{
			GzColor* colorPointer = (GzColor*)(valueList[0]); 
			render->flatcolor[0] = (*colorPointer)[0];
			render->flatcolor[1] = (*colorPointer)[1];
			render->flatcolor[2] = (*colorPointer)[2];
		}
		else if (nameList[i] == GZ_TEXTURE_MAP)
		{
			GzTexture _Tex = (GzTexture)(valueList[i]); 
			render->tex_fun = _Tex;
		}
		else if (nameList[i] == GZ_AASHIFTX)
		{
			float* shiftx = (float*)(valueList[i]);
			render->a_offsetX = (*shiftx);
		}
		else if (nameList[i] == GZ_AASHIFTY)
		{
			float* shifty = (float*)(valueList[i]);
			render->a_offsetY = (*shifty);
		}
	}
	return GZ_SUCCESS;
}

void ColorCalculate(GzRender *render, Normals _NormalTemp, GzColor Color)
{
	Vector3 _SpecularSum, _DiffuseSum;
	Vector3 _ColorOfLight;
	Vector3 _DirectionOfLight;
	Vector3 _NormalVector = Vector3(_NormalTemp.x, _NormalTemp.y, _NormalTemp.z);
	Vector3 _Edge = Vector3(0, 0, -1);
	 
	for(int i = 0; i < render->numlights; i++)
	{
		//CALCULATING THE COLOR OF LIGHT
		_ColorOfLight = Vector3(render->lights[i].color[0], render->lights[i].color[1], render->lights[i].color[2]);

		_DirectionOfLight = Vector3(render->lights[i].direction[0], render->lights[i].direction[1], render->lights[i].direction[2]);
		_DirectionOfLight.Normalise();
		Vector3 R = Vector3(0.0,0.0,0.0);

		if(_NormalVector.Dot(_DirectionOfLight) > 0 && _NormalVector.Dot(_Edge) < 0 || _NormalVector.Dot(_DirectionOfLight) < 0 && _NormalVector.Dot(_Edge) > 0 )
			continue;

		else if(_NormalVector.Dot(_DirectionOfLight) < 0 && _NormalVector.Dot(_Edge) < 0)
		{	
			_NormalVector= Vector3(-_NormalTemp.x, -_NormalTemp.y, -_NormalTemp.z);
			R = Vector3((2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_x - _DirectionOfLight.m_x), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_y - _DirectionOfLight.m_y), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_z - _DirectionOfLight.m_z));
		}
		else
		{
			R = Vector3((2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_x - _DirectionOfLight.m_x), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_y - _DirectionOfLight.m_y), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_z - _DirectionOfLight.m_z));
		}

		R.Normalise();
		float RdotE = R.Dot(_Edge);
		if(RdotE < 0)
			RdotE = 0;
		if(RdotE > 1)
			RdotE = 1.0;

		_SpecularSum.m_x +=  render->Ks[0] * _ColorOfLight.m_x * (pow(RdotE,render->spec));
		_SpecularSum.m_y += render->Ks[1] * _ColorOfLight.m_y * (pow(RdotE,render->spec));
		_SpecularSum.m_z += render->Ks[2] * _ColorOfLight.m_z * (pow(RdotE,render->spec));

		float NdotL = _NormalVector.Dot(_DirectionOfLight);
		_DiffuseSum.m_x += _ColorOfLight.m_x * NdotL;
		_DiffuseSum.m_y += _ColorOfLight.m_y * NdotL;
		_DiffuseSum.m_z += _ColorOfLight.m_z * NdotL;


	}
	Color[0] = _SpecularSum.m_x + render->Kd[0] * _DiffuseSum.m_x + render->Ka[0] * render->ambientlight.color[0];
	Color[1] = _SpecularSum.m_y + render->Kd[1] * _DiffuseSum.m_y + render->Ka[1] * render->ambientlight.color[1];
	Color[2] = _SpecularSum.m_z + render->Kd[2] * _DiffuseSum.m_z + render->Ka[2] * render->ambientlight.color[2];
	
	//Clamp colors

	if(Color[0] > 1.0)
		Color[0] =  1.0;

	if(Color[1] > 1.0)
		Color[1] =  1.0;

	if(Color[2] > 1.0)
		Color[2] =  1.0;

	if(Color[0] < 0.0)
		Color[0] =  0.0;

	if(Color[1] < 0.0)
		Color[1] =  0.0;

	if(Color[2] < 0.0)
		Color[2] =  0.0;
}

void ColorCalculateGoroud(GzRender *render, Normals _NormalTemp, GzColor Color)
{
	Vector3 _SpecularSum, _DiffuseSum;
	Vector3 _ColorOfLight;
	Vector3 _DirectionOfLight;
	Vector3 _NormalVector = Vector3(_NormalTemp.x, _NormalTemp.y, _NormalTemp.z);
	Vector3 _Edge = Vector3(0, 0, -1);
	 
	for(int i = 0; i < render->numlights; i++)
	{
		//CALCULATING THE COLOR OF LIGHT
		_ColorOfLight = Vector3(render->lights[i].color[0], render->lights[i].color[1], render->lights[i].color[2]);

		_DirectionOfLight = Vector3(render->lights[i].direction[0], render->lights[i].direction[1], render->lights[i].direction[2]);
		_DirectionOfLight.Normalise();
		Vector3 R = Vector3(0.0,0.0,0.0);

		if(_NormalVector.Dot(_DirectionOfLight) > 0 && _NormalVector.Dot(_Edge) < 0 || _NormalVector.Dot(_DirectionOfLight) < 0 && _NormalVector.Dot(_Edge) > 0 )
			continue;

		else if(_NormalVector.Dot(_DirectionOfLight) < 0 && _NormalVector.Dot(_Edge) < 0)
		{	
			_NormalVector= Vector3(-_NormalTemp.x, -_NormalTemp.y, -_NormalTemp.z);
			R = Vector3((2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_x - _DirectionOfLight.m_x), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_y - _DirectionOfLight.m_y), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_z - _DirectionOfLight.m_z));
		}
		else
		{
			R = Vector3((2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_x - _DirectionOfLight.m_x), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_y - _DirectionOfLight.m_y), (2*(_NormalVector.Dot(_DirectionOfLight))*_NormalVector.m_z - _DirectionOfLight.m_z));
		}

		R.Normalise();
		float RdotE = R.Dot(_Edge);
		if(RdotE < 0)
			RdotE = 0;
		if(RdotE > 1)
			RdotE = 1.0;

		_SpecularSum.m_x += _ColorOfLight.m_x * (pow(RdotE,render->spec));
		_SpecularSum.m_y += _ColorOfLight.m_y * (pow(RdotE,render->spec));
		_SpecularSum.m_z += _ColorOfLight.m_z * (pow(RdotE,render->spec));

		float NdotL = _NormalVector.Dot(_DirectionOfLight);
		_DiffuseSum.m_x += _ColorOfLight.m_x * NdotL;
		_DiffuseSum.m_y += _ColorOfLight.m_y * NdotL;
		_DiffuseSum.m_z += _ColorOfLight.m_z * NdotL;


	}
	Color[0] = _SpecularSum.m_x + _DiffuseSum.m_x + render->ambientlight.color[0];
	Color[1] = _SpecularSum.m_y + _DiffuseSum.m_y + render->ambientlight.color[1];
	Color[2] = _SpecularSum.m_z + _DiffuseSum.m_z + render->ambientlight.color[2];
	
	//Clamp colors

	if(Color[0] > 1.0)
		Color[0] =  1.0;

	if(Color[1] > 1.0)
		Color[1] =  1.0;

	if(Color[2] > 1.0)
		Color[2] =  1.0;

	if(Color[0] < 0.0)
		Color[0] =  0.0;

	if(Color[1] < 0.0)
		Color[1] =  0.0;

	if(Color[2] < 0.0)
		Color[2] =  0.0;
}

void ConvertToPerspective(Point *_Point)
{
	float zNew = _Point->z / INT_MAX - _Point->z;

	_Point->x = _Point->x/(zNew + 1);
	_Point->y = _Point->y/(zNew + 1);
	//_Point->z = _Point->z/(zNew + 1);
}
void ConvertToAffine(Point *_Point)
{
	float zNew = _Point->z;

	_Point->x = _Point->x * (zNew + 1);
	_Point->y = _Point->y * (zNew + 1);
	//_Point->z = _Point->z * (zNew + 1);
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
 	if (render == NULL || nameList == NULL || valueList == NULL)
		return GZ_FAILURE;

	if(render->interp_mode == GZ_FLAT)
	{
		
	}

	GzCoord* _Vert;
	GzCoord* _Normal;
	GzCoord* _NormalFlat = (GzCoord*)valueList[1];
	GzTextureIndex* _TexIndex;

	for (int i = 0; i < numParts; ++i)
	{
		if (nameList[i] == GZ_POSITION) 
			_Vert = (GzCoord*) valueList[i];

		else if (nameList[i] == GZ_NORMAL) 
			_Normal = (GzCoord*) valueList[i];

		else if(nameList[i] == GZ_TEXTURE_INDEX)
			_TexIndex = (GzTextureIndex*) valueList[i];
	}

	Point _VertN[3];
	for(int k=0 ; k<3 ; k++)
	{
		_VertN[k].x= _Vert[k][X];
		_VertN[k].y= _Vert[k][Y];
		_VertN[k].z= _Vert[k][Z];
	}

	Plane _FlatPlane;
	_FlatPlane.CalculateValues(_VertN);
	_FlatPlane.CalculateAndNormalize();
	
	Vector3 _FlatNormal = _FlatPlane.m_Normal;

	float a = _FlatNormal.m_x * render->Xnorm[normalMatrixlevel][0][0] + _FlatNormal.m_y * render->Xnorm[normalMatrixlevel][0][1] + _FlatNormal.m_z * render->Xnorm[normalMatrixlevel][0][2] + render->Xnorm[normalMatrixlevel][0][3];
	float b = _FlatNormal.m_x * render->Xnorm[normalMatrixlevel][1][0] + _FlatNormal.m_y * render->Xnorm[normalMatrixlevel][1][1] + _FlatNormal.m_z * render->Xnorm[normalMatrixlevel][1][2] + render->Xnorm[normalMatrixlevel][1][3];
	float c = _FlatNormal.m_x * render->Xnorm[normalMatrixlevel][2][0] + _FlatNormal.m_y * render->Xnorm[normalMatrixlevel][2][1] + _FlatNormal.m_z * render->Xnorm[normalMatrixlevel][2][2] + render->Xnorm[normalMatrixlevel][2][3];

	Normals _FlatShadingNormal;
	_FlatShadingNormal.x = a;
	_FlatShadingNormal.y = b;
	_FlatShadingNormal.z = c;

	GzCoord* _TransformTriangle = new GzCoord[3];

	GzMatrix _TriangleMatrix;
	GzMatrix _NormalMatrix;
	GzPixel _TempPixel;

	GzMatrix _WorldMatrix;

	_WorldMatrix[0][0] = _Vert[0][X];	
	_WorldMatrix[0][1] = _Vert[1][X];	
	_WorldMatrix[0][2] = _Vert[2][X];	
	_WorldMatrix[0][3] = 1;
	_WorldMatrix[1][0] = _Vert[0][Y];	
	_WorldMatrix[1][1] = _Vert[1][Y];	
	_WorldMatrix[1][2] = _Vert[2][Y];	
	_WorldMatrix[1][3] = 1;
	_WorldMatrix[2][0] = _Vert[0][Z];	
	_WorldMatrix[2][1] = _Vert[1][Z];	
	_WorldMatrix[2][2] = _Vert[2][Z];	
	_WorldMatrix[2][3] = 1;
	_WorldMatrix[3][0] = 1;				
	_WorldMatrix[3][1] = 1;				
	_WorldMatrix[3][2] = 1;				
	_WorldMatrix[3][3] = 1;

	GzMatrixMutiply(_TriangleMatrix,render->Ximage[render->matlevel],_WorldMatrix);

	_TransformTriangle[0][X] = _TriangleMatrix[0][0] / _TriangleMatrix[3][0];  
	_TransformTriangle[0][Y] = _TriangleMatrix[1][0] / _TriangleMatrix[3][0];
	_TransformTriangle[0][Z] = _TriangleMatrix[2][0] / _TriangleMatrix[3][0];

	_TransformTriangle[1][X] = _TriangleMatrix[0][1] / _TriangleMatrix[3][1];
	_TransformTriangle[1][Y] = _TriangleMatrix[1][1] / _TriangleMatrix[3][1];
	_TransformTriangle[1][Z] = _TriangleMatrix[2][1] / _TriangleMatrix[3][1];

	_TransformTriangle[2][X] = _TriangleMatrix[0][2] / _TriangleMatrix[3][2];
	_TransformTriangle[2][Y] = _TriangleMatrix[1][2] / _TriangleMatrix[3][2];
	_TransformTriangle[2][Z] = _TriangleMatrix[2][2] / _TriangleMatrix[3][2];

	for (int k = 0; k < 3; k++) 
	{ 
		_TransformTriangle[k][X] = _TransformTriangle[k][X] - render->a_offsetX ;  
		_TransformTriangle[k][Y] = _TransformTriangle[k][Y] - render->a_offsetY;
	}

	GzMatrix _WorldMatrixNormal;
	_WorldMatrixNormal[0][0] = _Normal[0][X];	
	_WorldMatrixNormal[0][1] = _Normal[1][X];	
	_WorldMatrixNormal[0][2] = _Normal[2][X];	
	_WorldMatrixNormal[0][3] = 1;

	_WorldMatrixNormal[1][0] = _Normal[0][Y];	
	_WorldMatrixNormal[1][1] = _Normal[1][Y];	
	_WorldMatrixNormal[1][2] = _Normal[2][Y];	
	_WorldMatrixNormal[1][3] = 1;

	_WorldMatrixNormal[2][0] = _Normal[0][Z];	
	_WorldMatrixNormal[2][1] = _Normal[1][Z];	
	_WorldMatrixNormal[2][2] = _Normal[2][Z];	
	_WorldMatrixNormal[2][3] = 1;

	_WorldMatrixNormal[3][0] =	1;				
	_WorldMatrixNormal[3][1] =	1;				
	_WorldMatrixNormal[3][2] =	1;				
	_WorldMatrixNormal[3][3] = 1;

	GzMatrixMutiply(_NormalMatrix,render->Xnorm[normalMatrixlevel],_WorldMatrixNormal);


	GzCoord* _Normalvertex = new GzCoord[3];
	_Normalvertex[0][X] = _NormalMatrix[0][0] / _NormalMatrix[3][0];  
	_Normalvertex[0][Y] = _NormalMatrix[1][0] / _NormalMatrix[3][0];
	_Normalvertex[0][Z] = _NormalMatrix[2][0] / _NormalMatrix[3][0];
					  				  
	_Normalvertex[1][X] = _NormalMatrix[0][1] / _NormalMatrix[3][1];
	_Normalvertex[1][Y] = _NormalMatrix[1][1] / _NormalMatrix[3][1];
	_Normalvertex[1][Z] = _NormalMatrix[2][1] / _NormalMatrix[3][1];
					  				  
	_Normalvertex[2][X] = _NormalMatrix[0][2] / _NormalMatrix[3][2];
	_Normalvertex[2][Y] = _NormalMatrix[1][2] / _NormalMatrix[3][2];
	_Normalvertex[2][Z] = _NormalMatrix[2][2] / _NormalMatrix[3][2];


	VectorNormalise(_Normalvertex[0]);
	VectorNormalise(_Normalvertex[1]);
	VectorNormalise(_Normalvertex[2]);




	//RASYERIZER
	Point _VertexPosition[3];									 
	Normals _Normal1[3];

	_VertexPosition[0].SetValues(_TransformTriangle[0][X], _TransformTriangle[0][Y], _TransformTriangle[0][Z]);
	_VertexPosition[1].SetValues(_TransformTriangle[1][X], _TransformTriangle[1][Y], _TransformTriangle[1][Z]);
	_VertexPosition[2].SetValues(_TransformTriangle[2][X], _TransformTriangle[2][Y], _TransformTriangle[2][Z]);													 			   

	_Normal1[0].SetValues(_Normalvertex[0][X], _Normalvertex[0][Y], _Normalvertex[0][Z]);
	_Normal1[1].SetValues(_Normalvertex[1][X], _Normalvertex[1][Y], _Normalvertex[1][Z]);
	_Normal1[2].SetValues(_Normalvertex[2][X], _Normalvertex[2][Y], _Normalvertex[2][Z]);

	//texture stuff
	Point _TextureCoords[3];
	for(int k=0 ; k<3 ; k++)
	{
		_TextureCoords[k].x= _TexIndex[k][X];
		_TextureCoords[k].y= _TexIndex[k][Y];
		_TextureCoords[k].z= _TransformTriangle[k][Z];
	}

	for(int i=0 ; i<3 ; i++)
	{
		float vz = _TextureCoords[i].z / (INT_MAX - _TextureCoords[i].z);
		_TextureCoords[i].x = _TextureCoords[i].x/(vz+1);
		_TextureCoords[i].y = _TextureCoords[i].y/(vz+1);
	}

	Point _TempPoint;
	Normals _TempNormal;
	Point _TempTexture;
	for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
      {
        if(_VertexPosition[i].y<_VertexPosition[j].y)
        {
          _TempPoint = _VertexPosition[i];
          _VertexPosition[i] = _VertexPosition[j];
          _VertexPosition[j] = _TempPoint;

		  _TempNormal = _Normal1[i];
          _Normal1[i] = _Normal1[j];
          _Normal1[j] = _TempNormal;

		  _TempPoint = _TextureCoords[i];
          _TextureCoords[i] = _TextureCoords[j];
          _TextureCoords[j] = _TempPoint;

        }      
      }      
    }

	

	GzColor colorvertex[3];
	if(render->interp_mode == GZ_COLOR)
	{
		ColorCalculateGoroud(render, _Normal1[0], colorvertex[0]);
		ColorCalculateGoroud(render, _Normal1[1], colorvertex[1]);
		ColorCalculateGoroud(render, _Normal1[2], colorvertex[2]);
	}
	else
	{
		ColorCalculate(render, _Normal1[0], colorvertex[0]);
		ColorCalculate(render, _Normal1[1], colorvertex[1]);
		ColorCalculate(render, _Normal1[2], colorvertex[2]);
	}
	

	float A1 = (_VertexPosition[1].y-_VertexPosition[0].y);
	float B1 = - (_VertexPosition[1].x-_VertexPosition[0].x);
	float C1 =  (_VertexPosition[1].x-_VertexPosition[0].x) * _VertexPosition[0].y  - (_VertexPosition[1].y-_VertexPosition[0].y) * _VertexPosition[0].x;

	float A2 = (_VertexPosition[2].y-_VertexPosition[1].y);
	float B2 = - (_VertexPosition[2].x-_VertexPosition[1].x);
	float C2 =  (_VertexPosition[2].x-_VertexPosition[1].x) * _VertexPosition[1].y  - (_VertexPosition[2].y-_VertexPosition[1].y) * _VertexPosition[1].x;//dXY - dYX

	float A3 = (_VertexPosition[0].y-_VertexPosition[2].y);
	float B3 = - (_VertexPosition[0].x-_VertexPosition[2].x);
	float C3 =  (_VertexPosition[0].x-_VertexPosition[2].x) * _VertexPosition[2].y  - (_VertexPosition[0].y-_VertexPosition[2].y) * _VertexPosition[2].x;//dXY - dYX

	int _RightEdge = -1; 

	float xP = (- C3 - B3 * _VertexPosition[1].y) / (A3);
	if(xP < _VertexPosition[1].x)
		_RightEdge = 1;
	else if (xP > _VertexPosition[1].x)
		_RightEdge = 3;

	float ymin = _VertexPosition[0].y;
	float ymax = _VertexPosition[2].y;
	float xmin = 10000;
	float xmax = 0;

	for(int i =0;i<3;i++)
	{
		if(xmin > _VertexPosition[i].x)
			xmin = _VertexPosition[i].x;

		if(xmax < _VertexPosition[i].x)
			xmax = _VertexPosition[i].x;
	}


	float E1,E2,E3 ; 	
	
	//PLane equation for Z 
	Plane _ZPlane;
	_ZPlane.CalculateValues(_VertexPosition);
	
	//PLane equation for Red 
	Plane _RedPlane;
	_RedPlane.CalculateColor(_VertexPosition, colorvertex, 0);

	//PLane equation for Green 
	Plane _GreenPlane;
	_GreenPlane.CalculateColor(_VertexPosition, colorvertex, 1);
	
	//PLane equation for Blue 
	Plane _BluePlane;
	_BluePlane.CalculateColor(_VertexPosition, colorvertex, 2);

	//PLane equation for Nx 
	Plane _NormalXPlane;
	_NormalXPlane.CalculateNormals(_VertexPosition, _Normal1, 0);

	//PLane equation for Ny 
	Plane _NormalYPlane;
	_NormalYPlane.CalculateNormals(_VertexPosition, _Normal1, 1);

	//PLane equation for Nz 
	Plane _NormalZPlane;
	_NormalZPlane.CalculateNormals(_VertexPosition, _Normal1, 2);

	//PLane equation for u
	Plane _TextureUPlane;
	_TextureUPlane.CalculateTextures(_VertexPosition, _TextureCoords, 0);

	//PLane equation for v
	Plane _TextureVPlane;
	_TextureVPlane.CalculateTextures(_VertexPosition, _TextureCoords, 1);

	for(int i=ceil(xmin);i<ceil(xmax);i++)
	{
		for(int j=ceil(ymin);j<ceil(ymax);j++)
		{
			E1 = (_VertexPosition[1].y-_VertexPosition[0].y) * (i - _VertexPosition[0].x) - (_VertexPosition[1].x - _VertexPosition[0].x) * (j - _VertexPosition[0].y);
			E2 = (_VertexPosition[2].y-_VertexPosition[1].y) * (i - _VertexPosition[1].x) - (_VertexPosition[2].x - _VertexPosition[1].x) * (j - _VertexPosition[1].y);
			E3 = (_VertexPosition[0].y-_VertexPosition[2].y) * (i - _VertexPosition[2].x) - (_VertexPosition[0].x - _VertexPosition[2].x) * (j - _VertexPosition[2].y);

			bool _Condition1 = ((E1 > 0 && E2 > 0  && E3 > 0) && (_RightEdge == 3));
			bool _Condition2 = ((E1 < 0 && E2 < 0  && E3 < 0) && (_RightEdge == 1));
			if(_Condition1 ||_Condition2)
			{
				float zval,rval,gval,bval;
				zval = (- _ZPlane.m_A * i - _ZPlane.m_B * j - _ZPlane.m_D) / _ZPlane.m_C ;

				GzTextureIndex _UV;
				_UV[0] = (- _TextureUPlane.m_A * i - _TextureUPlane.m_B * j - _TextureUPlane.m_D) / _TextureUPlane.m_C ;
				_UV[1] = (- _TextureVPlane.m_A * i - _TextureVPlane.m_B * j - _TextureVPlane.m_D) / _TextureVPlane.m_C ;


				GzTextureIndex _UVNew;
				float _NewZ = zval /(INT_MAX - zval);
				_UVNew[0] = _UV[0] * (_NewZ + 1);
				_UVNew[1] = _UV[1] * (_NewZ + 1);

				GzColor texColor;
				if (render->tex_fun != NULL)
				{		 
					render->tex_fun(_UVNew[0], _UVNew[1], texColor); 
				}

				if(render->interp_mode == GZ_COLOR)
				{
					rval = (- _RedPlane.m_A * i - _RedPlane.m_B * j - _RedPlane.m_D) / _RedPlane.m_C ;
					gval = (- _GreenPlane.m_A * i - _GreenPlane.m_B * j - _GreenPlane.m_D) / _GreenPlane.m_C ;
					bval = (- _BluePlane.m_A * i - _BluePlane.m_B * j - _BluePlane.m_D) / _BluePlane.m_C ; 

					rval *= texColor[0];
					gval *= texColor[1];
					bval *= texColor[2];

					if(rval > 1)
						rval = 1;
					if(gval > 1)
						gval = 1;
					if(bval > 1)
						bval = 1;
				}

				else if(render->interp_mode == GZ_NORMALS)
				{
					Normals PixelNormal;
					GzCoord PixelColor;

					PixelNormal.x = (- _NormalXPlane.m_A * i - _NormalXPlane.m_B * j - _NormalXPlane.m_D) / _NormalXPlane.m_C ;
					PixelNormal.y = (- _NormalYPlane.m_A * i - _NormalYPlane.m_B * j - _NormalYPlane.m_D) / _NormalYPlane.m_C ;
					PixelNormal.z = (- _NormalZPlane.m_A * i - _NormalZPlane.m_B * j - _NormalZPlane.m_D) / _NormalZPlane.m_C ;
					if (render->tex_fun != NULL)
				{	
					render->Ka[RED] = render->Kd[RED] = texColor[RED];
					render->Ka[GREEN] = render->Kd[GREEN] = texColor[GREEN];
					render->Ka[BLUE] = render->Kd[BLUE] = texColor[BLUE];
					}
					ColorCalculate(render, PixelNormal, PixelColor);
					if(PixelColor[0] > 1)
						PixelColor[0] = 1;
					if(PixelColor[1] > 1)
						PixelColor[1] = 1;
					if(PixelColor[2] >1)
						PixelColor[2] = 1;
					rval = PixelColor[0];
					gval = PixelColor[1];
					bval = PixelColor[2];

				}
				else if(render->interp_mode == GZ_FLAT)
				{	
					GzColor color;
					ColorCalculate(render, _FlatShadingNormal, color);
					rval = color[0];
					gval = color[1];
					bval = color[2];
				}
				GzGetDisplay(render->display, i, j, &_TempPixel.red, &_TempPixel.green, &_TempPixel.blue,&_TempPixel.alpha,&_TempPixel.z);
				
				if(zval < _TempPixel.z)
				{
					GzPutDisplay(render->display,i,j,ctoi(rval),ctoi(gval),ctoi(bval),1,zval);
				}
				

			}
			else if(E1 == 0 && (_RightEdge == 1))
			{
				float zval,rval,gval,bval;
				zval = (- _ZPlane.m_A * i - _ZPlane.m_B * j - _ZPlane.m_D) / _ZPlane.m_C ;

				GzTextureIndex _UV;
				_UV[0] = (- _TextureUPlane.m_A * i - _TextureUPlane.m_B * j - _TextureUPlane.m_D) / _TextureUPlane.m_C ;
				_UV[1] = (- _TextureVPlane.m_A * i - _TextureVPlane.m_B * j - _TextureVPlane.m_D) / _TextureVPlane.m_C ;


				GzTextureIndex _UVNew;
				float _NewZ = zval /(INT_MAX - zval);
				_UVNew[0] = _UV[0] * (_NewZ + 1);
				_UVNew[1] = _UV[1] * (_NewZ + 1);

				GzColor texColor;
				if (render->tex_fun != NULL)
				{		 
					render->tex_fun(_UVNew[0], _UVNew[1], texColor); 
				}

				if(render->interp_mode == GZ_COLOR)
				{
					rval = (- _RedPlane.m_A * i - _RedPlane.m_B * j - _RedPlane.m_D) / _RedPlane.m_C ;
					gval = (- _GreenPlane.m_A * i - _GreenPlane.m_B * j - _GreenPlane.m_D) / _GreenPlane.m_C ;
					bval = (- _BluePlane.m_A * i - _BluePlane.m_B * j - _BluePlane.m_D) / _BluePlane.m_C ; 

					rval *= texColor[0];
					gval *= texColor[1];
					bval *= texColor[2];

					if(rval > 1)
						rval = 1;
					if(gval > 1)
						gval = 1;
					if(bval > 1)
						bval = 1;
				}

				else if(render->interp_mode == GZ_NORMALS)
				{
					Normals PixelNormal;
					GzCoord PixelColor;

					PixelNormal.x = (- _NormalXPlane.m_A * i - _NormalXPlane.m_B * j - _NormalXPlane.m_D) / _NormalXPlane.m_C ;
					PixelNormal.y = (- _NormalYPlane.m_A * i - _NormalYPlane.m_B * j - _NormalYPlane.m_D) / _NormalYPlane.m_C ;
					PixelNormal.z = (- _NormalZPlane.m_A * i - _NormalZPlane.m_B * j - _NormalZPlane.m_D) / _NormalZPlane.m_C ;

					render->Ka[RED] = render->Kd[RED] = texColor[RED];
					render->Ka[GREEN] = render->Kd[GREEN] = texColor[GREEN];
					render->Ka[BLUE] = render->Kd[BLUE] = texColor[BLUE];

					ColorCalculate(render, PixelNormal, PixelColor);
					if(PixelColor[0] > 1)
						PixelColor[0] = 1;
					if(PixelColor[1] > 1)
						PixelColor[1] = 1;
					if(PixelColor[2] >1)
						PixelColor[2] = 1;
					rval = PixelColor[0];
					gval = PixelColor[1];
					bval = PixelColor[2];

				}
				else if(render->interp_mode == GZ_FLAT)
				{	GzColor color;
					ColorCalculate(render, _FlatShadingNormal, color);
					rval = color[0];
					gval = color[1];
					bval = color[2];
				}
				GzGetDisplay(render->display, i, j, &_TempPixel.red, &_TempPixel.green, &_TempPixel.blue,&_TempPixel.alpha,&_TempPixel.z);
				
				if(zval < _TempPixel.z)
					GzPutDisplay(render->display,i,j,ctoi(rval),ctoi(gval),ctoi(bval),1,zval);

			}
			else if(E3 == 0 && (_RightEdge == 3) )
			{
				float zval,rval,gval,bval;
				zval = (- _ZPlane.m_A * i - _ZPlane.m_B * j - _ZPlane.m_D) / _ZPlane.m_C ;

				GzTextureIndex _UV;
				_UV[0] = (- _TextureUPlane.m_A * i - _TextureUPlane.m_B * j - _TextureUPlane.m_D) / _TextureUPlane.m_C ;
				_UV[1] = (- _TextureVPlane.m_A * i - _TextureVPlane.m_B * j - _TextureVPlane.m_D) / _TextureVPlane.m_C ;


				GzTextureIndex _UVNew;
				float _NewZ = zval /(INT_MAX - zval);
				_UVNew[0] = _UV[0] * (_NewZ + 1);
				_UVNew[1] = _UV[1] * (_NewZ + 1);

				GzColor texColor;
				if (render->tex_fun != NULL)
				{		 
					render->tex_fun(_UVNew[0], _UVNew[1], texColor); 
				}

				if(render->interp_mode == GZ_COLOR)
				{
					rval = (- _RedPlane.m_A * i - _RedPlane.m_B * j - _RedPlane.m_D) / _RedPlane.m_C ;
					gval = (- _GreenPlane.m_A * i - _GreenPlane.m_B * j - _GreenPlane.m_D) / _GreenPlane.m_C ;
					bval = (- _BluePlane.m_A * i - _BluePlane.m_B * j - _BluePlane.m_D) / _BluePlane.m_C ; 

					rval *= texColor[0];
					gval *= texColor[1];
					bval *= texColor[2];

					if(rval > 1)
						rval = 1;
					if(gval > 1)
						gval = 1;
					if(bval > 1)
						bval = 1;
				}
				else if(render->interp_mode == GZ_NORMALS)
				{
					Normals PixelNormal;
					GzCoord PixelColor;

					PixelNormal.x = (- _NormalXPlane.m_A * i - _NormalXPlane.m_B * j - _NormalXPlane.m_D) / _NormalXPlane.m_C ;
					PixelNormal.y = (- _NormalYPlane.m_A * i - _NormalYPlane.m_B * j - _NormalYPlane.m_D) / _NormalYPlane.m_C ;
					PixelNormal.z = (- _NormalZPlane.m_A * i - _NormalZPlane.m_B * j - _NormalZPlane.m_D) / _NormalZPlane.m_C ;

					render->Ka[RED] = render->Kd[RED] = texColor[RED];
					render->Ka[GREEN] = render->Kd[GREEN] = texColor[GREEN];
					render->Ka[BLUE] = render->Kd[BLUE] = texColor[BLUE];

					ColorCalculate(render, PixelNormal, PixelColor);
					if(PixelColor[0] > 1)
						PixelColor[0] = 1;
					if(PixelColor[1] > 1)
						PixelColor[1] = 1;
					if(PixelColor[2] >1)
						PixelColor[2] = 1;
					rval = PixelColor[0];
					gval = PixelColor[1];
					bval = PixelColor[2];

				}
				else if(render->interp_mode == GZ_FLAT)
				{	
					GzColor color;
					ColorCalculate(render, _FlatShadingNormal, color);
					rval = color[0];
					gval = color[1];
					bval = color[2];
				}
				GzGetDisplay(render->display, i, j, &_TempPixel.red, &_TempPixel.green, &_TempPixel.blue,&_TempPixel.alpha,&_TempPixel.z);
				
				if(zval < _TempPixel.z)
				{
					GzPutDisplay(render->display,i,j,ctoi(rval),ctoi(gval),ctoi(bval),1,zval);
				}
			}
		}
	}
	return GZ_SUCCESS;
}



