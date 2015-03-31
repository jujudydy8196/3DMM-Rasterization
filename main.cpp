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
    "/Users/judy/Documents/senior/3DMM/3dmm_hw1/model/quad.obj"
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
		cout << "vs_z: " << vs_z << "ve_z: " << ve_z << endl;
		cout << "connect v" << vStart << "(" << vs_x << "," << vs_y << "," << vs_z
				<< ")-v" << vEnd << "(" << ve_x << "," << ve_y << "," << ve_z << endl;
		if (ve_x == vs_x)
				for (int y = min(vs_y,ve_y)+1; y<max(vs_y,ve_y); y++)
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
						cout << "z: " << z << endl;
						framebuffer.draw(vs_x,y,z,vec3(z/300.f));
				}
		else
		{
				m = (float(ve_y-vs_y))/(float(ve_x-vs_x));
				b = vs_y - (m*vs_x);
				//cout << "m: " << m << " b: " << b << endl;	
				if (abs(m)<1)
				{
					 cout << "m<1" << endl;
					 int x = min(vs_x,ve_x)+1;
					 float y;
					 for ( x; x<max(vs_x,ve_x); x++)
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
							 cout << "z: " << z << endl;
							 //cout << "draw: (" << x << "," << int(y) << ")" << endl; 
							 //cout << 1.f/f\loat(vs_z) << endl;
							 framebuffer.draw(x,int(y),z,vec3(z/300.f));
					 }
				}
				else
				{
						cout << "m>1" << endl;
						int y = min(vs_y,ve_y)+1;
						float x;
						//cout << "y init: " << y << endl;
						for (y; y<max(ve_y,vs_y); y++)
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
							  cout << "z: " << z << endl;
							//cout << "draw: (" << int(x) << "," << y << ")" << endl; 
								//cout << 1.f/float(vs_z) << endl;
								framebuffer.draw(int(x),y,z,vec3(z/300.f));
						}
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
		// draw the framebuffer
		ix = int(((1.f-margin)*curX)*screenSide) + screenWidth_half - 1;
		iy = int(((1.f-margin)*curY)*screenSide) + screenHeight_half - 1;
		iz = int(((1.f-margin)*curZ)*screenSide) + max(screenWidth,screenHeight)/2 - 1;
		modelPtr[curModelIdx]->projects[3*i] = ix;
		modelPtr[curModelIdx]->projects[3*i+1] = iy;
		modelPtr[curModelIdx]->projects[3*i+2] = iz;
		cout << "ix: " << ix << " iy: " << iy << " iz: " << iz <<" iz/300: " << iz/300.f <<  endl;
	  if (iz/300.f)
		framebuffer.draw(ix, iy, curZ, vec3(iz/300));
		}

	for(int j=0; j<modelPtr[curModelIdx]->numTriangles; j++)
	{
		v0 = modelPtr[curModelIdx]->triangles[j].vIndices[0];
		v1 = modelPtr[curModelIdx]->triangles[j].vIndices[1]; 
		v2 = modelPtr[curModelIdx]->triangles[j].vIndices[2]; 
		//cout << "triangles index: " << v0 << " " << v1 << " " << v2 << endl;
	
		DrawLine(v0,v1);
		DrawLine(v1,v2);
		DrawLine(v2,v0);

	}


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
