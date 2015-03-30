#include "framebuffer.h"

Framebuffer::Framebuffer(int w, int h)
  : width(w), height(h), clearColor(0.f)
{
	int numPixel = w*h;
	colorBuffer = new vec3[numPixel];
	depthBuffer = new float[numPixel];
	clear();
}

Framebuffer::~Framebuffer() 
{
	if(colorBuffer) delete[] colorBuffer;
	if(depthBuffer) delete[] depthBuffer;
}

void Framebuffer::clear() const 
{
	int numPixel = width * height;
	for (int i = 0; i < numPixel; ++i) {
		colorBuffer[i] = clearColor;
		depthBuffer[i] = DEPTH_INF;
	}
}

void Framebuffer::setClearColor(const vec3 color) 
{
	clearColor = color;
}

void Framebuffer::draw(int ix, int iy, float z, vec3 c) const {
	if(ix < 0 || ix >= width || iy < 0 || iy >= height) {
		return;
	}
	int idx = iy * width + ix;
	if (z <= depthBuffer[idx]) {
		colorBuffer[idx] = c;
		depthBuffer[idx] = z;
	}
}

const vec3* Framebuffer::getPixels() const {
	return colorBuffer;
}

void Framebuffer::writePPM(string filename) const {

	char* buffer = new char[height*width*3];
	for(int i = 0; i < height*width; ++i) {
		buffer[3*i  ] = (char)(colorBuffer[i].x*255.f);
		buffer[3*i+1] = (char)(colorBuffer[i].y*255.f);
		buffer[3*i+2] = (char)(colorBuffer[i].z*255.f);
	}

	ofstream ofs(filename.c_str(), ios::binary);
	ofs << "P6" << endl
		<< width << " " << height << endl
		<< "255" << endl;
	ofs.write(buffer, height*width*3);
	ofs.close();

	delete [] buffer;
}
