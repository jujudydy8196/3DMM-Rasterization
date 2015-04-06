#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#endif

//#include <GL/glut.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include "model.h"
#include "framebuffer.h"
using namespace std;

#define PI 3.14159265
#define ROUND(x) (int)(x+0.5)
Model** modelPtr;

/* frame size */
int screenWidth = 800;
int screenHeight = 600;
int screenWidth_half = 400;
int screenHeight_half = 300;

/* theta */
double rotateSpeed = 0.05;
double thetaX = 0;
double thetaY = 0;

int curModelIdx;
bool culling;
Framebuffer framebuffer(screenWidth, screenHeight);

/* model names */
const char* modelNames[] = {
    "/Users/judy/Documents/senior/3DMM/3dmm_hw1/model/couch.obj"
//	"model/quad.obj"
//	"model/couch.obj",
//	"model/blaze.obj",
//	"model/ateneal.obj",
//	"model/venusm.obj",
//	"model/bunnyC.obj",
//	"model/duck4KN.obj",
//	"model/happy10KN.obj",
//	"model/dragon10KN.obj",
//	"model/elephant16KN.obj",
//	"model/Statue_of_Liberty.obj",
//	"model/Nissan_Pathfinder.obj"
};
const int numModels = sizeof(modelNames) / sizeof(char*);

void printHelp() 
{
	printf("===========================================\n");
	printf("  H/h: Show help menu                      \n");
	printf("  M/m: Select model                        \n");
	printf("  UP/DOWN: Rotate along x-axis             \n");
	printf("  LEFT/RIGHT: Rotate along y-axis          \n");
	printf("  C/c: Toggle back-face culling            \n");
	printf("  B/b: Toggle background color             \n");
	printf("  S/s: Save image                          \n");
	printf("  Q/q: Quit                                \n");
	printf("===========================================\n\n");
}

void init() 
{
	modelPtr = new Model*[numModels];
	for(int i = 0; i < numModels; ++i) {
		/* load models */
		modelPtr[i] = readObj(modelNames[i]);

		/* normalize coordinates to [-0.5, +0.5] */
		float xmax, xmin, ymax, ymin, zmax, zmin, longest;
		float xmean, ymean, zmean;
		xmax = xmin = modelPtr[i]->vertices[3];
		ymax = ymin = modelPtr[i]->vertices[4];
		zmax = zmin = modelPtr[i]->vertices[5];
		for(int v = 2; v <= modelPtr[i]->numVertices; ++v) {
			xmax = max(xmax, modelPtr[i]->vertices[3*v]);
			xmin = min(xmin, modelPtr[i]->vertices[3*v]);
			ymax = max(ymax, modelPtr[i]->vertices[3*v+1]);
			ymin = min(ymin, modelPtr[i]->vertices[3*v+1]);
			zmax = max(zmax, modelPtr[i]->vertices[3*v+2]);
			zmin = min(zmin, modelPtr[i]->vertices[3*v+2]);
		}
		longest = max(xmax - xmin, ymax - ymin);
		longest = max(longest, zmax - zmin);
		xmean = 0.5f*(xmax + xmin);
		ymean = 0.5f*(ymax + ymin);
		zmean = 0.5f*(zmax + zmin);
		
		for(int v = 1; v <= modelPtr[i]->numVertices; ++v) {
			modelPtr[i]->vertices[3*v  ] = (modelPtr[i]->vertices[3*v  ] - xmean)/longest;
			modelPtr[i]->vertices[3*v+1] = (modelPtr[i]->vertices[3*v+1] - ymean)/longest;
			modelPtr[i]->vertices[3*v+2] = (modelPtr[i]->vertices[3*v+2] - zmean)/longest;
		}
	}

	/* initialize parameters */
	curModelIdx = 0;
	culling = true;
}

void DrawLine(int vStart, int vEnd)
{
		int vs_x,vs_y,vs_z,ve_x,ve_y,ve_z;
		float m,b,z;
		vs_x = modelPtr[curModelIdx]->projects[3*vStart  ];
		vs_y = modelPtr[curModelIdx]->projects[3*vStart+1];
		vs_z = modelPtr[curModelIdx]->projects[3*vStart+2];
		ve_x = modelPtr[curModelIdx]->projects[3*vEnd  ];
		ve_y = modelPtr[curModelIdx]->projects[3*vEnd+1];
		ve_z = modelPtr[curModelIdx]->projects[3*vEnd+2];
		int dis_x = abs(vs_x-ve_x);
		int dis_y = abs(vs_y-ve_y);
		int dis_z = ve_z-vs_z;
		//cout << "vs_z: " << vs_z << "ve_z: " << ve_z << endl;
		//cout << "connect v" << vStart << "(" << vs_x << "," << vs_y << "," << vs_z
				//<< ")-v" << vEnd << "(" << ve_x << "," << ve_y << "," << ve_z << endl;
		if (ve_x == vs_x)
				for (int y = min(vs_y,ve_y); y<=max(vs_y,ve_y); y++)
				{
						if (dis_z == 0)
								z = vs_z;
						else
						{
								if (min(vs_y,ve_y)==vs_y )
								{
										//cout << vs_z << " " <<  float(dis_z) << " " << float(y-min(vs_y,ve_y)) / dis_y << " " << float(dis_z)* float(y-min(vs_y,ve_y)) / dis_y << endl;
										z = vs_z + (float(dis_z) * float(y-min(vs_y,ve_y)) / dis_y );
								}
								else
								{
										//cout << ve_z << " " <<  float(dis_z) << " " << float(y-min(vs_y,ve_y)) / dis_y << " " << float(dis_z)* float(y-min(vs_y,ve_y)) / dis_y << endl;
										z = ve_z - (float(dis_z) * float(y-min(vs_y,ve_y)) / dis_y );
								}
						}
						//cout << "z: " << z << endl;
						framebuffer.draw(vs_x,y,z,vec3(1.f-(z/670.f)));
				}
		else
		{
				m = (float(ve_y-vs_y))/(float(ve_x-vs_x));
				b = vs_y - (m*vs_x);
				//cout << "m: " << m << " b: " << b << endl;	
				if (abs(m)<1)
				{
					 //cout << "m<1" << endl;
					 int x = min(vs_x,ve_x);//+1;
					 float y;
					 for ( x; x<=max(vs_x,ve_x); x++)
					 {
							 y = m*x+b;
							 if (dis_z == 0)
									 z = vs_z;
							 else
							 {
									 if (min(vs_x,ve_x) == vs_x)
											 z = vs_z + dis_z * ( float(x-min(vs_x,ve_x)) / dis_x );
									 else
											 z = ve_z -dis_z * ( float(x-min(vs_x,ve_x)) / dis_x );
							 }
							 //cout << "z: " << z << endl;
							 //cout << "draw: (" << x << "," << int(y) << ")" << endl; 
							 //cout << 1.f/f\loat(vs_z) << endl;
							 framebuffer.draw(x,int(y),z,vec3(1.f-(z/670.f)));
					 }
				}
				else
				{
						//cout << "m>1" << endl;
						int y = min(vs_y,ve_y);//+1;
						float x;
						//cout << "y init: " << y << endl;
						for (y; y<=max(ve_y,vs_y); y++)
						{
								//cout << "y: " << y << endl;
								x = (float(y-b))/m;
								if (dis_z ==0 )
										z = vs_z;
								else
								{
										if (min(vs_y,ve_y)==vs_y)
										{
												//cout << vs_z << dis_z << float(y-min(vs_y,ve_y)) / dis_y << endl;
												z = vs_z + (dis_z * ( float(y-min(vs_y,ve_y)) / dis_y ));
										}
										else
										{
												//cout << ve_z << dis_z << float(y-min(vs_y,ve_y)) / dis_y << endl;
												z = ve_z - (dis_z * ( float(y-min(vs_y,ve_y)) / dis_y ));
										}
								}
								//cout << "z: " << z << endl;
							//cout << "draw: (" << int(x) << "," << y << ")" << endl; 
								//cout << 1.f/float(vs_z) << endl;
								framebuffer.draw(int(x),y,z,vec3(1.f-(z/670.f)));
						}
				}
		}
}

void drawline(int vStart[3] , int vEnd[3])
{
		float m,b,z;
		int dis_x = abs(vStart[0]-vEnd[0]);
		int dis_y = abs(vStart[1]-vEnd[1]);
		int dis_z = vEnd[2]-vStart[2];
		//cout << "vs_z: " << vs_z << "ve_z: " << ve_z << endl;
		//cout << "connect v" <<  "(" << vStart[0] << "," << vStart[1] << "," << vStart[2] 
				//<< ")-v(" << vEnd[0] << "," << vEnd[1] << "," << vEnd[2] << endl;
		if (vEnd[0] == vStart[0])
				for (int y = min(vStart[1],vEnd[1]); y<=max(vStart[1],vEnd[1]); y++)
				{
						if (dis_z == 0)
								z = vStart[2];
						else
						{
								if (min(vStart[1],vEnd[1])==vStart[1] )
								{
										//cout << vs_z << " " <<  float(dis_z) << " " << float(y-min(vs_y,ve_y)) / dis_y << " " << float(dis_z)* float(y-min(vs_y,ve_y)) / dis_y << endl;
										z = vStart[2] + (float(dis_z) * float(y-min(vStart[1],vEnd[1])) / dis_y );
								}
								else
								{
										//cout << ve_z << " " <<  float(dis_z) << " " << float(y-min(vs_y,ve_y)) / dis_y << " " << float(dis_z)* float(y-min(vs_y,ve_y)) / dis_y << endl;
										z = vEnd[2] - (float(dis_z) * float(y-min(vStart[1],vEnd[1])) / dis_y );
								}
						}
						//cout << "z: " << z << endl;
						framebuffer.draw(vStart[0],y,z,vec3(1.f-(z/670.f)));
				}
		else
		{
				m = (float(vEnd[1]-vStart[1]))/(float(vEnd[0]-vStart[0]));
				b = vStart[1] - (m*vStart[0]);
				//cout << "m: " << m << " b: " << b << endl;	
				if (abs(m)<1)
				{
					 //cout << "m<1" << endl;
					 int x = min(vStart[0],vEnd[0]);//+1;
					 float y;
					 for ( x; x<=max(vStart[0],vEnd[0]); x++)
					 {
							 y = m*x+b;
							 if (dis_z == 0)
									 z = vStart[2];
							 else
							 {
									 if (min(vStart[0],vEnd[0]) == vStart[0])
											 z = vStart[2] + dis_z * ( float(x-min(vStart[0],vEnd[0])) / dis_x );
									 else
											 z = vEnd[2] -dis_z * ( float(x-min(vStart[0],vEnd[0])) / dis_x );
							 }
							 //cout << "z: " << z << endl;
							 //cout << "draw: (" << x << "," << int(y) << ")" << endl; 
							 //cout << 1.f/f\loat(vs_z) << endl;
							 framebuffer.draw(x,ROUND(y),z,vec3(1.f-(z/670.f)));
					 }
				}
				else
				{
						//cout << "m>1" << endl;
						int y = min(vStart[1],vEnd[1]);//+1;
						float x;
						//cout << "y init: " << y << endl;
						for (y; y<=max(vStart[1],vEnd[1]); y++)
						{
								//cout << "y: " << y << endl;
								x = (float(y-b))/m;
								if (dis_z ==0 )
										z = vStart[2];
								else
								{
										if (min(vEnd[1],vStart[1])==vStart[1])
										{
												//cout << vs_z << dis_z << float(y-min(vs_y,ve_y)) / dis_y << endl;
												z = vStart[2] + (dis_z * ( float(y-min(vStart[1],vEnd[1])) / dis_y ));
										}
										else
										{
												//cout << ve_z << dis_z << float(y-min(vs_y,ve_y)) / dis_y << endl;
												z = vEnd[2] - (dis_z * ( float(y-min(vStart[1],vEnd[1])) / dis_y ));
										}
								}
								//cout << "z: " << z << endl;
							//cout << "draw: (" << int(x) << "," << y << ")" << endl; 
								//cout << 1.f/float(vs_z) << endl;
								framebuffer.draw(ROUND(x),y,z,vec3(1.f-(z/670.f)));
						}
				}
		}
}
void fillTopTriangle ( int sharp[3], int v1[3], int v2[3])
{
		//cout << " top " << sharp[0] << " " << sharp[1] << " " << v1[0] << " " << v1[1] << " " << v2[0] << " " << v2[1] <<  endl;
		//float ms1 = float( sharp[0] - v1[0] ) / float( sharp[1]-v1[1]);
		//float ms2 = float( sharp[0] - v2[0] ) / float( sharp[1]-v2[1]);
		//int nowx1 = sharp[0];
		//int nowx2 = sharp[0];

		float m1,m2,b1,b2;
		m1 = float(sharp[1]-v1[1])/float(sharp[0]-v1[0]);
		m2 = float(sharp[1]-v2[1])/float(sharp[0]-v2[0]);
		b1 = float(v1[1]) - m1*(float)v1[0];
		b2 = float(v2[1]) - m2*(float)v2[0];

		//cout << "m1: " << m1 << " m2: " << m2 << " b1: " << b1 << " b2: " << b2 << endl;

		//cout << "ms1: " << float(sharp[0]-v1[0]) << "/" << float(sharp[1]-v1[1]) <<" " << ms1
				//<< " ms2: " << float(sharp[0]-v2[0]) << "/" << float(sharp[1]-v2[1]) <<" " << ms2 << endl;

		for (int y = sharp[1]; y <= v1[1] ; y++ )
		{
			  int x1 = (sharp[0]-v1[0] == 0 ) ? sharp[0] : ROUND((float)(y-b1)/m1);
			  int x2 = (sharp[0]-v2[0] == 0 ) ? sharp[0] : ROUND((float)(y-b2)/m2);
				int L1[3] = {x1,y,sharp[2]+ ROUND( (v1[2]-sharp[2]) * (y-sharp[1])/(v1[1]-sharp[1]))};
				int L2[3] = {x2,y,sharp[2]+ ROUND( (v2[2]-sharp[2]) * (y-sharp[1])/(v2[1]-sharp[1]))};
				//cout << "L1: (" << L1[0] << "," << L1[1] << "," << L1[2] << ")" << endl;
				//cout << "L2: (" << L2[0] << "," << L2[1] << "," << L2[2] << ")" << endl;
				drawline(L1,L2);
	
				//nowx1 = (roundf)(nowx1+ms1);
				//nowx2 = (roundf)(nowx2+ms2);
		}
}
void fillBottomTriangle( int sharp[3], int v1[3], int v2[3])
{
		//cout << " Botton! " << endl;
		//float ms1 = float( sharp[0] - v1[0] ) / float( sharp[1]-v1[1]);
		//float ms2 = float( sharp[0] - v2[0] ) / float( sharp[1]-v2[1]);
		//int nowx1 = sharp[0];
		//int nowx2 = sharp[0];

		//cout << "(" << sharp[0] << "," << sharp[1] << ")   (" << v1[0] << "," << v1[1] << ")   ("<< v2[0] << "," << v2[1] << ")" << endl;

		//cout << "ms1: " << ms1 << " ms2: " << ms2 << endl;
		float m1,m2,b1,b2;
		m1 = float(sharp[1]-v1[1])/float(sharp[0]-v1[0]);
		m2 = float(sharp[1]-v2[1])/float(sharp[0]-v2[0]);
		b1 = float(v1[1]) - m1*(float)v1[0];
		b2 = float(v2[1]) - m2*(float)v2[0];
		//cout << "m1: " << m1 << " m2: " << m2 << " b1: " << b1 << " b2: " << b2 << endl;


		for (int y = sharp[1]; y >= v1[1] ; y-- )
		{
			  int x1 = (sharp[0]-v1[0] == 0 ) ? sharp[0] : ROUND((float)(y-b1)/m1);
			  int x2 = (sharp[0]-v2[0] == 0 ) ? sharp[0] : ROUND((float)(y-b2)/m2);

				int L1[3] = {x1,y,sharp[2]+ ROUND( (v1[2]-sharp[2]) * (y-sharp[1])/(v1[1]-sharp[1]))};
				int L2[3] = {x2,y,sharp[2]+ ROUND( (v2[2]-sharp[2]) * (y-sharp[1])/(v2[1]-sharp[1]))};
				//cout << "L1: (" << L1[0] << "," << L1[1] << "," << L1[2] << ")" << endl;
				//cout << "L2: (" << L2[0] << "," << L2[1] << "," << L2[2] << ")" << endl;

				drawline(L1,L2);
				//nowx1 = (roundf)(nowx1-ms1);
				//nowx2 = (roundf)(nowx2-ms2);
		}

}
void fillTriangles ( Triangle* t )
{
		int vv0,vv1,vv2;
		vv0 = t->vIndices[0];
		vv1 = t->vIndices[1];
		vv2 = t->vIndices[2];

		//cout << "vv0: " << vv0 << " vv1: " << vv1 << " vv2: " << vv2 << endl;

		int v[3][3] = {
				{ modelPtr[curModelIdx]->projects[3*vv0], modelPtr[curModelIdx]->projects[3*vv0+1],modelPtr[curModelIdx]->projects[3*vv0+2]},
				{ modelPtr[curModelIdx]->projects[3*vv1], modelPtr[curModelIdx]->projects[3*vv1+1],modelPtr[curModelIdx]->projects[3*vv1+2]},
				{ modelPtr[curModelIdx]->projects[3*vv2], modelPtr[curModelIdx]->projects[3*vv2+1],modelPtr[curModelIdx]->projects[3*vv2+2]} };

		//float vec1[3] = {float(v[1][0]-v[0][0]),float(v[1][1]-v[0][1]),float(v[1][2]-v[0][2])};
		//float vec2[3] = {float(v[2][0]-v[0][0]),float(v[2][1]-v[0][1]),float(v[2][2]-v[0][2])};

		//float mag1 = sqrt((vec1[0]*vec1[0])+(vec1[1]*vec1[1])+(vec1[2]*vec1[2]));
		//float mag2 = sqrt((vec2[0]*vec2[0])+(vec2[1]*vec2[1])+(vec2[2]*vec2[2]));
	
		////cout << "mag1: " << mag1 << " mag2: " << mag2 << endl;
		//for (int i=0 ; i<3; i++)
		//{
				////cout << "bvec1: " << vec1[i] << " bvec2: " << vec2[i] << endl;
				//vec1[i] /= mag1;
				//vec2[i] /= mag2;
				////cout << "vec1: " << vec1[i] << " vec2: " << vec2[i] << endl;
		//}
		
		////for back culling
		//t->faceNormal[0] = (vec1[1]*vec2[2])-(vec1[2]*vec2[1]);
		//t->faceNormal[1] = (vec1[2]*vec2[0])-(vec1[0]*vec2[2]);
		//t->faceNormal[2] = (vec1[0]*vec2[1])-(vec1[1]*vec2[0]);

		//float dotProduct = float(v[0][0])*(t->faceNormal[0]) + float(v[0][1])*(t->faceNormal[1]) + float(v[0][2]-1)*(t->faceNormal[2]);
		////cout << "dotProduct: " << dotProduct << endl;
		////cout << "---------------------" << endl;
		//if ( dotProduct < 0.0f )
		//{
				////cout << "no display" << endl;
				//return;
		//}
		int ymax = max(v[0][1],max(v[1][1],v[2][1]));
		int ysorted[3];
		if (ymax == v[0][1])
		{
				ysorted[0] =0;
				ysorted[1] = (v[1][1] >= v[2][1])? 1 : 2;
				ysorted[2] = (v[1][1] < v[2][1])? 1 : 2;
		}
		else if (ymax == v[1][1])
		{
				ysorted[0] =1;
				ysorted[1] = (v[0][1] >= v[2][1])? 0 : 2;
				ysorted[2] = (v[0][1] < v[2][1])? 0 : 2;
		}
		else
		{
				ysorted[0] =2;
				ysorted[1] = (v[1][1] >= v[0][1])? 1:0;
				ysorted[2] = (v[1][1] < v[0][1])? 1:0;
		}

		//cout << v[0][1] << " " << v[1][1] << " " << v[2][1] << endl;
		//cout << "ysorted index: " ;
		//for (int i=0; i<3;i++)
				//cout << ysorted[i] << " " ;
		//cout << endl;
				DrawLine(vv0,vv1);
				DrawLine(vv1,vv2);
				DrawLine(vv2,vv0);
	
		if (v[ysorted[0]][1] != v[ysorted[2]][1] )
		{	
				if (v[ysorted[1]][1] == v[ysorted[2]][1] )
						fillBottomTriangle(v[ysorted[0]],v[ysorted[1]],v[ysorted[2]]);
				else if (v[ysorted[0]][1] == v[ysorted[1]][1])
						fillTopTriangle(v[ysorted[2]],v[ysorted[0]],v[ysorted[1]]);
				else
				{
						int v4[3] = { ROUND(v[ysorted[0]][0] + ((float)(v[ysorted[1]][1] - v[ysorted[0]][1]) / (float)( v[ysorted[2]][1] - v[ysorted[0]][1] )) * ( v[ysorted[2]][0] - v[ysorted[0]][0])) , v[ysorted[1]][1] , ROUND(v[ysorted[0]][2] + ((float)(v[ysorted[1]][1] - v[ysorted[0]][1]) / (float)( v[ysorted[2]][1] - v[ysorted[0]][1] )) * ( v[ysorted[2]][2] - v[ysorted[0]][2]))};
						//cout << "v4: " << v4[0] << " " << v4[1] << " " << v4[2] << endl; 
						fillBottomTriangle(v[ysorted[0]],v[ysorted[1]],v4);
						fillTopTriangle(v[ysorted[2]],v[ysorted[1]],v4);
				}
		}

}

void displayFunc() 
{
	framebuffer.clear();

	/* draw point cloud */
	/* TODO: change the following code to do rasterization */
	int ix, iy, iz, v0, v1, v2;
	int screenSide = min(screenWidth, screenHeight);
	float margin = 0.1f;
	float prevX, prevY, prevZ, curX, curY, curZ;
	//float v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, m01, m12, m20, b;
	for(int i = 1; i <= modelPtr[curModelIdx]->numVertices; ++i) {
		prevX = modelPtr[curModelIdx]->vertices[3*i  ];
		prevY = modelPtr[curModelIdx]->vertices[3*i+1];
		prevZ = modelPtr[curModelIdx]->vertices[3*i+2];
		//cout << "prev " << i << " : " << prevX << " " << prevY << " " << prevZ << endl;
		// rotate along x-axis
		curY = prevY*cos(thetaX)-prevZ*sin(thetaX);
		curZ = prevY*sin(thetaX)+prevZ*cos(thetaX);
		// prevY = curY;
		prevZ = curZ;
		// rotate along y-axis
		curX = prevZ*sin(thetaY)+prevX*cos(thetaY);
		curZ = prevZ*cos(thetaY)-prevX*sin(thetaY);
		modelPtr[curModelIdx]->curs[3*i] = curX;
		modelPtr[curModelIdx]->curs[3*i+1] = curY;
		modelPtr[curModelIdx]->curs[3*i+2] = curZ;

		// draw the framebuffer
		ix = int(((1.f-margin)*curX)*screenSide) + screenWidth_half - 1;
		iy = int(((1.f-margin)*curY)*screenSide) + screenHeight_half - 1;
		iz = int(((1.f-margin)*curZ)*screenSide) + max(screenWidth,screenHeight)/2 - 1;
		modelPtr[curModelIdx]->projects[3*i] = ix;
		modelPtr[curModelIdx]->projects[3*i+1] = iy;
		modelPtr[curModelIdx]->projects[3*i+2] = iz;
		//cout << "ix: " << ix << " iy: " << iy << " iz: " << iz <<" iz/300: " << iz/670.f <<  endl;
		//if (iz/670.f)
			//framebuffer.draw(ix, iy, curZ, vec3(1.f-(iz/670.f)));
		}

	for(int j=0; j<modelPtr[curModelIdx]->numTriangles; j++)
	{
		v0 = modelPtr[curModelIdx]->triangles[j].vIndices[0];
		v1 = modelPtr[curModelIdx]->triangles[j].vIndices[1]; 
		v2 = modelPtr[curModelIdx]->triangles[j].vIndices[2]; 
		//cout << "triangles index: " << v0 << " " << v1 << " " << v2 << endl;
	
		float vt[3][3] = {
				{ modelPtr[curModelIdx]->curs[3*v0], modelPtr[curModelIdx]->curs[3*v0+1],modelPtr[curModelIdx]->curs[3*v0+2]},
				{ modelPtr[curModelIdx]->curs[3*v1], modelPtr[curModelIdx]->curs[3*v1+1],modelPtr[curModelIdx]->curs[3*v1+2]},
				{ modelPtr[curModelIdx]->curs[3*v2], modelPtr[curModelIdx]->curs[3*v2+1],modelPtr[curModelIdx]->curs[3*v2+2]} };
		float vec1[3] = {vt[1][0]-vt[0][0],vt[1][1]-vt[0][1],vt[1][2]-vt[0][2]};
		float vec2[3] = {vt[2][0]-vt[0][0],vt[2][1]-vt[0][1],vt[2][2]-vt[0][2]};

		float mag1 = sqrt((vec1[0]*vec1[0])+(vec1[1]*vec1[1])+(vec1[2]*vec1[2]));
		float mag2 = sqrt((vec2[0]*vec2[0])+(vec2[1]*vec2[1])+(vec2[2]*vec2[2]));
	
		//cout << "mag1: " << mag1 << " mag2: " << mag2 << endl;
		for (int i=0 ; i<3; i++)
		{
				//cout << "bvec1: " << vec1[i] << " bvec2: " << vec2[i] << endl;
				vec1[i] /= mag1;
				vec2[i] /= mag2;
				//cout << "vec1: " << vec1[i] << " vec2: " << vec2[i] << endl;
		}
		
		//Triangle* t = modelPtr[curModelIdx]->triangles[j];
		//for back culling
		modelPtr[curModelIdx]->triangles[j].faceNormal[0] = (vec1[1]*vec2[2])-(vec1[2]*vec2[1]);
		modelPtr[curModelIdx]->triangles[j].faceNormal[1] = (vec1[2]*vec2[0])-(vec1[0]*vec2[2]);
		modelPtr[curModelIdx]->triangles[j].faceNormal[2] = (vec1[0]*vec2[1])-(vec1[1]*vec2[0]);

		float dotProduct = float(vt[0][0])*(modelPtr[curModelIdx]->triangles[j].faceNormal[0]) + float(vt[0][1])*(modelPtr[curModelIdx]->triangles[j].faceNormal[1]) + float(vt[0][2]-3)*(modelPtr[curModelIdx]->triangles[j].faceNormal[2]);
		//cout << "dotProduct: " << dotProduct << endl;
		//cout << "---------------------" << endl;
		if (!culling)
			fillTriangles(&modelPtr[curModelIdx]->triangles[j]);
		else if ( dotProduct > 0.0f && culling)
		{
				//cout << "no display" << endl;
			fillTriangles(&modelPtr[curModelIdx]->triangles[j]);
		}

	}

	//for(int j=0; j<modelPtr[curModelIdx]->numTriangles;j++)
			//fillTriangles(&modelPtr[curModelIdx]->triangles[j]);


    /* display */
	glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, (const GLvoid*)framebuffer.getPixels());
	glutSwapBuffers();

	/* FPS counter */
	static clock_t prev = clock();
	static clock_t curr;
	static char title[20];
	static double refreshTime = 0.5;
	static int count = 0;
	++count;
	curr = clock();
	double t = (double)(curr - prev)/(double)CLOCKS_PER_SEC;
	if (t > refreshTime) {
		sprintf(title, "3DMM HW#1 Rasterization: %.2lf fps",  (double)count/t);
		glutSetWindowTitle(title);
		prev = curr;
		count = 0;
	}
}


void idleFunc() 
{
	glutPostRedisplay();
}

void keyboardFunc(unsigned char key, int x, int y) 
{
	switch (key) {
        
	// Quit
	case 'q':
	case 'Q':
		exit(0);
		break;
	// Help
	case 'h':
	case 'H':
		printHelp();
		break;
	// Save image
	case 's':
	case 'S':
		static time_t t;
		static char name[80];
		time(&t);
		strftime(name, sizeof(name), "%Y%m%d%H%M%S.ppm", localtime(&t));
		printf("Save framebuffer to %s\n", name);
		framebuffer.writePPM(name);
		break;
	// Select model
	case 'm':
		curModelIdx = (curModelIdx == numModels - 1)? 0 : (curModelIdx + 1);
		printf("Switch to model \"%s\"\n", modelNames[curModelIdx]);
		break;
	case 'M':
		curModelIdx = (curModelIdx == 0)? (numModels - 1) : (curModelIdx - 1);
		printf("Switch to model \"%s\"\n", modelNames[curModelIdx]);
		break;
	// Back-face culling
	case 'c':
	case 'C':
		culling = !culling;
		break;
	// Background color
	case 'b':
	case 'B':
		static bool isBlack = true;
		isBlack = !isBlack;
		framebuffer.setClearColor(isBlack? vec3(0.f) : vec3(1.f));
		break;
	// You can add more functions as you like.
	}
	glutPostRedisplay();
}

// Special keys
void specialFunc(int key, int x, int y) 
{
	switch (key) {
	case GLUT_KEY_RIGHT:
		thetaY = (thetaY > PI*2)? 0 : (thetaY + rotateSpeed);
		break;
	case GLUT_KEY_LEFT:
		thetaY = (thetaY < 0)? (PI*2) : (thetaY - rotateSpeed);
		break;
	case GLUT_KEY_UP:
		thetaX = (thetaX < 0)? (PI*2) : (thetaX - rotateSpeed);
		break;
	case GLUT_KEY_DOWN:
		thetaX = (thetaX > PI*2)? 0 : (thetaX + rotateSpeed);
		break;
	case GLUT_KEY_PAGE_UP:
    
		break;
	case GLUT_KEY_PAGE_DOWN:
    
		break;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(screenWidth, screenHeight);
	glutCreateWindow("3DMM HW#1 Rasterization");

	init();

	glutDisplayFunc(displayFunc);
	glutIdleFunc(idleFunc);
	// glutMouseFunc(mouseFunc);
	glutKeyboardFunc(keyboardFunc);
	glutSpecialFunc(specialFunc);

	glutMainLoop();
	
	for(int i = 0; i < numModels; ++i) {
		if(modelPtr[i]) delete [] modelPtr[i];
	}
	if(modelPtr) delete [] modelPtr;

	return 0;
}
