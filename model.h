#ifndef _MODEL_H_
#define _MODEL_H_

#include <iostream>
using namespace std;

class Triangle {
public:
	Triangle() {
		faceNormal = new float[3];
		vIndices = new int[3];
		nIndices = new int[3];
	}
	~Triangle() {
		delete [] vIndices;
		delete [] nIndices;
	}
	int* vIndices;	/* array of triangle vertex indices */
	int* nIndices;	/* array of triangle normal indices */
	float* faceNormal;
};

class Model {
public:
	Model(): numVertices(0), numNormals(0), numTriangles(0), vertices(nullptr), normals(nullptr), projects(nullptr), curs(nullptr), triangles(nullptr) {}
	~Model() { 
		if(vertices) delete [] vertices;
		if(normals) delete [] normals;
		if(triangles) delete [] triangles;
		if(projects) delete [] projects;
		if(curs) delete [] curs;
	}
	/* data members are declared public for easy access */
	int numVertices;		/* number of vertices in model */
	int numNormals;			/* number of normals in model */
	int numTriangles;		/* number of triangles in model */
	float* vertices;		/* array of vertices */
	float* normals;			/* array of normals */
	int*   projects;
	float* curs;
	Triangle* triangles;	/* array of triangles */
};

Model* readObj(const string filename);

#endif // _MODEL_H_
