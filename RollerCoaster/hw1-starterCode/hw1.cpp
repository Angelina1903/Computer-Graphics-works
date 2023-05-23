/*
  CSCI 420 Computer Graphics, USC
  Assignment 1: Height Fields with Shaders.
  C++ starter code

  Student username: liangang
*/

#include "basicPipelineProgram.h"
#include "openGLMatrix.h"
#include "imageIO.h"
#include "openGLHeader.h"
#include "glutHeader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openGLHeader.h"
#include "imageIO.h"
#include <iostream>
#include <cstring>
#include <vector>

#if defined(WIN32) || defined(_WIN32)
#ifdef _DEBUG
#pragma comment(lib, "glew32d.lib")
#else
#pragma comment(lib, "glew32.lib")
#endif
#endif

//#if defined(WIN32) || defined(_WIN32)
//char shaderBasePath[1024] = SHADER_BASE_PATH;
//#else
//char shaderBasePath[1024] = "../openGLHelper-starterCode";
//#endif

char shaderBasePath[1024] = "../openGLHelper-starterCode";
char textureBasePath[1024] = "../openGLHelper-texture";

using namespace std;

int mousePos[2]; // x,y coordinate of the mouse position

int leftMouseButton = 0; // 1 if pressed, 0 if not 
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0; // 1 if pressed, 0 if not

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;

// state of the world
float landRotate[3] = { 0.0f, 0.0f, 0.0f };
float landTranslate[3] = { 0.0f, 0.0f, 0.0f };
float landScale[3] = { 1.0f, 1.0f, 1.0f };

int windowWidth = 1280;
int windowHeight = 720;
char windowTitle[512] = "CSCI 420 homework II";

ImageIO* heightmapImage;

float *lineVertices, *trackColors;
GLuint triVertexBuffer, triColorVertexBuffer;
GLuint triVertexArray;
int sizeTri;

OpenGLMatrix matrix;
BasicPipelineProgram* pipelineProgram;


GLint h_modelViewMatrix, h_projectionMatrix;

// for spline calculation
float step = 0.001f;
float s = 0.5f;
glm::mat4 basis(
    glm::vec4(-s, 2.0f*s, -s, 0.0f),
    glm::vec4(2.0f-s, s-3.0f, 0.0f, 1.0f),
    glm::vec4(s-2.0f, 3.0f-2.0f*s, s, 0.0f),
    glm::vec4(s, -s, 0.0f, 0.0f));

// for line data
vector<glm::vec3> line_points;
vector<glm::vec3> line_tans;
vector<glm::vec3> line_binormals;
vector<glm::vec3> line_normals;

// for camera movements
vector<glm::vec3> cam_points;
vector<glm::vec3> cam_tans;
vector<glm::vec3> cam_binormals;
vector<glm::vec3> cam_normals;
int cam_step = 0;
float h_max = -100000000.0f;

// test
bool stopMove = false;
bool acc = false;
int screenShotCount = 0;

// texture handling
BasicPipelineProgram* textureProgram;
GLuint texHandle;
float *groundVertices, *groundTexVertices;
GLuint groundVertexBuffer, groundTexVertexBuffer;
GLuint groundVertexArray;
GLint ground_modelViewMatrix, ground_projectionMatrix;

// represents one control point along the spline 
struct Point
{
    double x;
    double y;
    double z;
};

// spline struct 
// contains how many control points the spline has, and an array of control points 
struct Spline
{
    int numControlPoints;
    Point* points;
};

// the spline array 
Spline* splines;
// total number of splines 
int numSplines;

int loadSplines(char* argv)
{
    char* cName = (char*)malloc(128 * sizeof(char));
    FILE* fileList;
    FILE* fileSpline;
    int iType, i = 0, j, iLength;

    // load the track file 
    fileList = fopen(argv, "r");
    if (fileList == NULL)
    {
        printf("can't open file\n");
        exit(1);
    }

    // stores the number of splines in a global variable 
    fscanf(fileList, "%d", &numSplines);

    splines = (Spline*)malloc(numSplines * sizeof(Spline));

    // reads through the spline files 
    for (j = 0; j < numSplines; j++)
    {
        i = 0;
        fscanf(fileList, "%s", cName);
        fileSpline = fopen(cName, "r");

        if (fileSpline == NULL)
        {
            printf("can't open file\n");
            exit(1);
        }

        // gets length for spline file
        fscanf(fileSpline, "%d %d", &iLength, &iType);

        // allocate memory for all the points
        splines[j].points = (Point*)malloc(iLength * sizeof(Point));
        splines[j].numControlPoints = iLength;

        // saves the data to the struct
        while (fscanf(fileSpline, "%lf %lf %lf",
            &splines[j].points[i].x,
            &splines[j].points[i].y,
            &splines[j].points[i].z) != EOF)
        {
            i++;
        }
    }

    free(cName);

    return 0;
}

int initTexture(const char* imageFilename, GLuint textureHandle)
{
    // read the texture image
    ImageIO img;
    ImageIO::fileFormatType imgFormat;
    ImageIO::errorType err = img.load(imageFilename, &imgFormat);

    if (err != ImageIO::OK)
    {
        printf("Loading texture from %s failed.\n", imageFilename);
        return -1;
    }

    // check that the number of bytes is a multiple of 4
    if (img.getWidth() * img.getBytesPerPixel() % 4)
    {
        printf("Error (%s): The width*numChannels in the loaded image must be a multiple of 4.\n", imageFilename);
        return -1;
    }

    // allocate space for an array of pixels
    int width = img.getWidth();
    int height = img.getHeight();
    unsigned char* pixelsRGBA = new unsigned char[4 * width * height]; // we will use 4 bytes per pixel, i.e., RGBA

    // fill the pixelsRGBA array with the image pixels
    memset(pixelsRGBA, 0, 4 * width * height); // set all bytes to 0
    for (int h = 0; h < height; h++)
        for (int w = 0; w < width; w++)
        {
            // assign some default byte values (for the case where img.getBytesPerPixel() < 4)
            pixelsRGBA[4 * (h * width + w) + 0] = 0; // red
            pixelsRGBA[4 * (h * width + w) + 1] = 0; // green
            pixelsRGBA[4 * (h * width + w) + 2] = 0; // blue
            pixelsRGBA[4 * (h * width + w) + 3] = 255; // alpha channel; fully opaque

            // set the RGBA channels, based on the loaded image
            int numChannels = img.getBytesPerPixel();
            for (int c = 0; c < numChannels; c++) // only set as many channels as are available in the loaded image; the rest get the default value
                pixelsRGBA[4 * (h * width + w) + c] = img.getPixel(w, h, c);
        }

    // bind the texture
    glBindTexture(GL_TEXTURE_2D, textureHandle);

    // initialize the texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelsRGBA);

    // generate the mipmaps for this texture
    glGenerateMipmap(GL_TEXTURE_2D);

    // set the texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // query support for anisotropic texture filtering
    GLfloat fLargest;
    glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
    printf("Max available anisotropic samples: %f\n", fLargest);
    // set anisotropic texture filtering
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 0.5f * fLargest);

    // query for any errors
    GLenum errCode = glGetError();
    if (errCode != 0)
    {
        printf("Texture initialization error. Error code: %d.\n", errCode);
        return -1;
    }

    // de-allocate the pixel array -- it is no longer needed
    delete[] pixelsRGBA;

    return 0;
}

// enable multitexturing
void setTextureUnit(GLint unit)
{
	glActiveTexture(unit); // select texture unit affected by subsequent texture calls
	// get a handle to the “textureImage” shader variable
	GLint h_textureImage = glGetUniformLocation(textureProgram->GetProgramHandle(), "textureImage");
	// deem the shader variable “textureImage” to read from texture unit “unit”
	glUniform1i(h_textureImage, unit - GL_TEXTURE0);
}

// write a screenshot to the specified filename
void saveScreenshot(const char* filename)
{
    unsigned char* screenshotData = new unsigned char[windowWidth * windowHeight * 3];
    glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

    ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

    if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
        cout << "File " << filename << " saved successfully." << endl;
    else cout << "Failed to save file " << filename << '.' << endl;

    delete[] screenshotData;
}

void displayFunc()
{
    // render some stuff...
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.LoadIdentity();

	/*glm::vec3 curr_pos = line_points[cam_step] + 0.2f*line_normals[cam_step];
	matrix.LookAt(curr_pos.x, curr_pos.y, curr_pos.z, 
		line_points[cam_step].x + line_tans[cam_step].x, line_points[cam_step].y + line_tans[cam_step].y, line_points[cam_step].z + line_tans[cam_step].z,
		line_normals[cam_step].x, line_normals[cam_step].y, line_normals[cam_step].z);*/

	glm::vec3 curr_pos = cam_points[cam_step] + 0.2f*cam_normals[cam_step];
	matrix.LookAt(curr_pos.x, curr_pos.y, curr_pos.z,
		cam_points[cam_step].x + cam_tans[cam_step].x, cam_points[cam_step].y + cam_tans[cam_step].y, cam_points[cam_step].z + cam_tans[cam_step].z,
		cam_normals[cam_step].x, cam_normals[cam_step].y, cam_normals[cam_step].z);

	//matrix.Translate(landTranslate[0], landTranslate[1], landTranslate[2]);
	//matrix.Scale(landScale[0], landScale[1], landScale[2]);
	//matrix.Rotate(landRotate[0], 1, 0, 0);
	//matrix.Rotate(landRotate[1], 0, 1, 0);
	//matrix.Rotate(landRotate[2], 0, 0, 1);

		// bind shader
	pipelineProgram->Bind();

	float view[16];
	matrix.GetMatrix(view);
	// get a handle to the program
	GLuint program1 = pipelineProgram->GetProgramHandle();
	// get a handle to the viewLightDirection shader variable
	GLint h_viewLightDirection =
		glGetUniformLocation(program1, "viewLightDirection");
	float lightDirection[3] = { 0, 1, 0 }; // the “Sun” at noon
	float viewLightDirection[3]; // light direction in the view space
	// the following line is pseudo-code:
	// viewLightDirection = (view * float4(lightDirection, 0.0)).xyz;
	viewLightDirection[0] = (view[0] * lightDirection[0]) + (view[4] * lightDirection[1]) + (view[8] * lightDirection[2]);
	viewLightDirection[1] = (view[1] * lightDirection[0]) + (view[5] * lightDirection[1]) + (view[9] * lightDirection[2]);
	viewLightDirection[2] = (view[2] * lightDirection[0]) + (view[6] * lightDirection[1]) + (view[10] * lightDirection[2]);
	// upload viewLightDirection to the GPU
	glUniform3fv(h_viewLightDirection, 1, viewLightDirection);

	float La[4] = { 1, 1, 1}; float Ld[4] = { 1, 1, 1}; float Ls[4] = { 1, 1, 1};
	float ka[4] = { 0.56, 0.56, 0.56}; float kd[4] = { 0.5, 0.5, 0.5}; float ks[4] = { 0, 0, 0};

	GLint h_La = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "La");
	GLint h_Ld = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "Ld");
	GLint h_Ls = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "Ls");
	glUniform4fv(h_La, 1, La); glUniform4fv(h_Ld, 1, Ld); glUniform4fv(h_Ls, 1, Ls);

	GLint h_ka = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "ka");
	GLint h_kd = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "kd");
	GLint h_ks = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "ks");
	GLint h_alpha = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "alpha");
	glUniform4fv(h_ka, 1, ka); glUniform4fv(h_kd, 1, kd); glUniform4fv(h_ks, 1, ks); glUniform1f(h_alpha, 1.0f);



	// continue with model transformations
	//openGLMatrix->Translate(x, y, z);

	// get a handle to the program
	GLuint program2 = pipelineProgram->GetProgramHandle();
	// get a handle to the normalMatrix shader variable
	GLint h_normalMatrix =
		glGetUniformLocation(program2, "normalMatrix");
	float n[16];
	matrix.SetMatrixMode(OpenGLMatrix::ModelView);

	matrix.GetNormalMatrix(n); // get normal matrix
	// upload n to the GPU
	GLboolean isRowMajor = GL_FALSE;
	glUniformMatrix4fv(h_normalMatrix, 1, isRowMajor, n);

    float m[16];
    matrix.SetMatrixMode(OpenGLMatrix::ModelView);
    matrix.GetMatrix(m);
    glUniformMatrix4fv(h_modelViewMatrix, 1, GL_FALSE, m);

    float p[16];
    matrix.SetMatrixMode(OpenGLMatrix::Projection);
    matrix.GetMatrix(p);
    glUniformMatrix4fv(h_projectionMatrix, 1, GL_FALSE, p);
    //
    //// bind shader
    //pipelineProgram->Bind();

    // set variable
    pipelineProgram->SetModelViewMatrix(m);
    pipelineProgram->SetProjectionMatrix(p);

    glBindVertexArray(triVertexArray);
    glDrawArrays(GL_TRIANGLES, 0, sizeTri*2);

	// texture setups
	glUniformMatrix4fv(ground_modelViewMatrix, 1, GL_FALSE, m);
	glUniformMatrix4fv(ground_projectionMatrix, 1, GL_FALSE, p);
	textureProgram->Bind();
	textureProgram->SetModelViewMatrix(m);
	textureProgram->SetProjectionMatrix(p);
	// select the active texture unit
	setTextureUnit(GL_TEXTURE0); // it is safe to always use GL_TEXTURE0
	// select the texture to use (“texHandle” was generated by glGenTextures)
	glBindTexture(GL_TEXTURE_2D, texHandle);
	// bind vao
	glBindVertexArray(groundVertexArray);
	glDrawArrays(GL_TRIANGLES, 0, 30);

    glutSwapBuffers();
}

void idleFunc()
{
    // do some stuff... 

    // for example, here, you can save the screenshots to disk (to make the animation)

	if (cam_step == cam_points.size() - 1) {
		cam_step = 0;
	}
	if (!stopMove) {
		if (acc && cam_step+3 < cam_points.size() - 1) cam_step += 3;
		else cam_step++;
	}// test

	//if (cam_step == line_points.size() - 1) {
	//	cam_step = 0;
	//}
	//if (!stopMove) {
	//	if (acc && cam_step + 5 < line_points.size() - 1) cam_step += 5;
	//	else cam_step++;
	//}// test

	//// screensots
	//if (screenShotCount <= 500) {
	//	string filename = "";
	//	if (screenShotCount < 10) {
	//		filename = "anim/00" + to_string(screenShotCount) + ".jpg";
	//	}
	//	else if (screenShotCount >= 10 && screenShotCount < 100) {
	//		filename = "anim/0" + to_string(screenShotCount) + ".jpg";
	//	}
	//	else if (screenShotCount >= 100) {
	//		filename = "anim/" + to_string(screenShotCount) + ".jpg";
	//	}
	//	saveScreenshot(filename.c_str());
	//	screenShotCount++;
	//}
    // make the screen update 
    glutPostRedisplay();
}

void reshapeFunc(int w, int h)
{
    glViewport(0, 0, w, h);

    matrix.SetMatrixMode(OpenGLMatrix::Projection);
    matrix.LoadIdentity();
    matrix.Perspective(54.0f, (float)w / (float)h, 0.01f, 1000.0f);
}

void mouseMotionDragFunc(int x, int y)
{
    // mouse has moved and one of the mouse buttons is pressed (dragging)

    // the change in mouse position since the last invocation of this function
    int mousePosDelta[2] = { x - mousePos[0], y - mousePos[1] };

	//if (!stopMove) return;//test

    switch (controlState)
    {
        // translate the landscape
    case TRANSLATE:
        if (leftMouseButton)
        {
            // control x,y translation via the left mouse button
            landTranslate[0] += mousePosDelta[0] * 0.01f;
            landTranslate[1] -= mousePosDelta[1] * 0.01f;
        }
        if (middleMouseButton)
        {
            // control z translation via the middle mouse button
            landTranslate[2] += mousePosDelta[1] * 0.01f;
        }
        break;

        // rotate the landscape
    case ROTATE:
        if (leftMouseButton)
        {
            // control x,y rotation via the left mouse button
            landRotate[0] += mousePosDelta[1];
            landRotate[1] += mousePosDelta[0];
        }
        if (middleMouseButton)
        {
            // control z rotation via the middle mouse button
            landRotate[2] += mousePosDelta[1];
        }
        break;

        // scale the landscape
    case SCALE:
        if (leftMouseButton)
        {
            // control x,y scaling via the left mouse button
            landScale[0] *= 1.0f + mousePosDelta[0] * 0.01f;
            landScale[1] *= 1.0f - mousePosDelta[1] * 0.01f;
        }
        if (middleMouseButton)
        {
            // control z scaling via the middle mouse button
            landScale[2] *= 1.0f - mousePosDelta[1] * 0.01f;
        }
        break;
    }

    // store the new mouse position
    mousePos[0] = x;
    mousePos[1] = y;
}

void mouseMotionFunc(int x, int y)
{
    // mouse has moved
    // store the new mouse position
    mousePos[0] = x;
    mousePos[1] = y;
}

void mouseButtonFunc(int button, int state, int x, int y)
{
    // a mouse button has has been pressed or depressed

    // keep track of the mouse button state, in leftMouseButton, middleMouseButton, rightMouseButton variables
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        leftMouseButton = (state == GLUT_DOWN);
        break;

    case GLUT_MIDDLE_BUTTON:
        middleMouseButton = (state == GLUT_DOWN);
        break;

    case GLUT_RIGHT_BUTTON:
        rightMouseButton = (state == GLUT_DOWN);
        break;
    }

    // keep track of whether CTRL and SHIFT keys are pressed
    switch (glutGetModifiers())
    {
    case GLUT_ACTIVE_CTRL:
        controlState = TRANSLATE;
        break;

    case GLUT_ACTIVE_SHIFT:
        controlState = SCALE;
        break;

        // if CTRL and SHIFT are not pressed, we are in rotate mode
    default:
        controlState = ROTATE;
        break;
    }

    // store the new mouse position
    mousePos[0] = x;
    mousePos[1] = y;
}

void keyboardFunc(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27: // ESC key
        exit(0); // exit the program
        break;

    case ' ':
        cout << "You pressed the spacebar." << endl;
		stopMove = !stopMove;//test
        break;

	case 'a':
		cout << "You pressed the ctrl key." << endl;
		acc = !acc;//test
		break;

    case 'x':
        // take a screenshot
        saveScreenshot("screenshot.jpg");
		//saveScreenshot("anim/testscreenshot.jpg");
        break;
    }
}

void initLineData() {
	// tangent vectors
	for (int i = 0; i < numSplines; ++i) {
		// for control points 1 to n-3
		for (int j = 1; j <= splines[i].numControlPoints - 3; ++j) {
			glm::mat3x4 control_mat(
				splines[i].points[j - 1].x, splines[i].points[j].x, splines[i].points[j + 1].x, splines[i].points[j + 2].x,
				splines[i].points[j - 1].y, splines[i].points[j].y, splines[i].points[j + 1].y, splines[i].points[j + 2].y,
				splines[i].points[j - 1].z, splines[i].points[j].z, splines[i].points[j + 1].z, splines[i].points[j + 2].z
			);
			for (float u = 0.0f; u <= 1.0f; u += step) {
				glm::vec4 u_vec(static_cast<float>(pow(u, 3)), static_cast<float>(pow(u, 2)), u, 1.0f);
				glm::vec3 p = u_vec * basis * control_mat;
				line_points.push_back(p);

				glm::vec4 tan_u_vec(3.0f*static_cast<float>(pow(u, 2)), 2.0f*u, 1.0f, 0.0f);
				glm::vec3 tan = glm::normalize(tan_u_vec * basis * control_mat);
				line_tans.push_back(tan);
			}
		}
	}

	// normal vectors & binormal vectors
	glm::vec3 v0(0.0f, 1.0f, 0.0f);

		// N0 = unit(T0 x V) and B0 = unit(T0 x N0)
	glm::vec3 n0 = glm::normalize(glm::cross(line_tans[0], v0));
	line_normals.push_back(n0);
	glm::vec3 b0 = glm::normalize(glm::cross(line_tans[0], n0));
	line_binormals.push_back(b0);

	for (size_t i = 1; i < line_tans.size(); ++i) {
		// N1 = unit(B0 x T1) and B1 = unit(T1 x N1)
		glm::vec3 ni = glm::normalize(glm::cross(line_binormals[i-1], line_tans[i]));
		line_normals.push_back(ni);
		glm::vec3 bi = glm::normalize(glm::cross(line_tans[i], ni));
		line_binormals.push_back(bi);
	}
}

void initCamData() {

	for (int i = 0; i < numSplines; ++i) {
		for (int j = 0; j < splines[i].numControlPoints; ++j) {
			if (splines[i].points[j].y > h_max) h_max = splines[i].points[j].y;
		}
	}
	// tangent vectors
	for (int i = 0; i < numSplines; ++i) {
		// for control points 1 to n-3
		for (int j = 1; j <= splines[i].numControlPoints - 3; ++j) {
			glm::mat3x4 control_mat(
				splines[i].points[j - 1].x, splines[i].points[j].x, splines[i].points[j + 1].x, splines[i].points[j + 2].x,
				splines[i].points[j - 1].y, splines[i].points[j].y, splines[i].points[j + 1].y, splines[i].points[j + 2].y,
				splines[i].points[j - 1].z, splines[i].points[j].z, splines[i].points[j + 1].z, splines[i].points[j + 2].z
			);
			float u = 0.0f;
			while (u <= 1.0f) {
				glm::vec4 u_vec(static_cast<float>(pow(u, 3)), static_cast<float>(pow(u, 2)), u, 1.0f);
				glm::vec3 p = u_vec * basis * control_mat;
				cam_points.push_back(p);

				glm::vec4 tan_u_vec(3.0f*static_cast<float>(pow(u, 2)), 2.0f*u, 1.0f, 0.0f);
				float dp_du = glm::length(tan_u_vec);

				glm::vec3 tan = glm::normalize(tan_u_vec * basis * control_mat);
				cam_tans.push_back(tan);

				//u += (0.0002f*sqrt(2.0f*9.8f*(abs(h_max - p.y))) / dp_du);
				u += step;
			}
		}
	}

	// normal vectors & binormal vectors
	glm::vec3 v0(0.0f, 1.0f, 0.0f);

	// N0 = unit(T0 x V) and B0 = unit(T0 x N0)
	glm::vec3 n0 = glm::normalize(glm::cross(cam_tans[0], v0));
	cam_normals.push_back(n0);
	glm::vec3 b0 = glm::normalize(glm::cross(cam_tans[0], n0));
	cam_binormals.push_back(b0);

	for (size_t i = 1; i < cam_tans.size(); ++i) {
		// N1 = unit(B0 x T1) and B1 = unit(T1 x N1)
		glm::vec3 ni = glm::normalize(glm::cross(cam_binormals[i - 1], cam_tans[i]));
		cam_normals.push_back(ni);
		glm::vec3 bi = glm::normalize(glm::cross(cam_tans[i], ni));
		cam_binormals.push_back(bi);
	}
}

void initRailCrossSectionData() {
	float alpha = 0.03f;
	int idx = 0;
	int idxx = 0;
	for (size_t i = 0; i < line_points.size() - 1; ++i) {
		// p, n, and b
		glm::vec3 p0 = line_points[i];
		glm::vec3 p1 = line_points[i + 1];

		glm::vec3 n0 = line_normals[i];
		glm::vec3 n1 = line_normals[i + 1];

		glm::vec3 b0 = line_binormals[i];
		glm::vec3 b1 = line_binormals[i + 1];

		// v0-v7
		glm::vec3 v0 = p0 + alpha * (-n0 + b0);
		glm::vec3 v1 = p0 + alpha * (n0 + b0);
		glm::vec3 v2 = p0 + alpha * (n0 - b0);
		glm::vec3 v3 = p0 + alpha * (-n0 - b0);

		glm::vec3 v4 = p1 + alpha * (-n1 + b1);
		glm::vec3 v5 = p1 + alpha * (n1 + b1);
		glm::vec3 v6 = p1 + alpha * (n1 - b1);
		glm::vec3 v7 = p1 + alpha * (-n1 - b1);

		// normals
		glm::vec3 n_up = n0;
		glm::vec3 n_down = -n0;
		glm::vec3 n_left = -b0;
		glm::vec3 n_right = b0;

		// setup vertices data array
		// v0,v1,v5 #1
		lineVertices[idx] = v0.x;
		lineVertices[idx + 1] = v0.y;
		lineVertices[idx + 2] = v0.z;
		idx += 3;
		lineVertices[idx] = v1.x;
		lineVertices[idx + 1] = v1.y;
		lineVertices[idx + 2] = v1.z;
		idx += 3;
		lineVertices[idx] = v5.x;
		lineVertices[idx + 1] = v5.y;
		lineVertices[idx + 2] = v5.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_right.x;
			trackColors[idxx + 1] = n_right.y;
			trackColors[idxx + 2] = n_right.z;
			idxx += 3;
		}
		// v0,v4,v5	#2
		lineVertices[idx] = v0.x;
		lineVertices[idx + 1] = v0.y;
		lineVertices[idx + 2] = v0.z;
		idx += 3;
		lineVertices[idx] = v4.x;
		lineVertices[idx + 1] = v4.y;
		lineVertices[idx + 2] = v4.z;
		idx += 3;
		lineVertices[idx] = v5.x;
		lineVertices[idx + 1] = v5.y;
		lineVertices[idx + 2] = v5.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_right.x;
			trackColors[idxx + 1] = n_right.y;
			trackColors[idxx + 2] = n_right.z;
			idxx += 3;
		}
		// v2,v1,v5	#3
		lineVertices[idx] = v2.x;
		lineVertices[idx + 1] = v2.y;
		lineVertices[idx + 2] = v2.z;
		idx += 3;
		lineVertices[idx] = v1.x;
		lineVertices[idx + 1] = v1.y;
		lineVertices[idx + 2] = v1.z;
		idx += 3;
		lineVertices[idx] = v5.x;
		lineVertices[idx + 1] = v5.y;
		lineVertices[idx + 2] = v5.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_up.x;
			trackColors[idxx + 1] = n_up.y;
			trackColors[idxx + 2] = n_up.z;
			idxx += 3;
		}
		// v2,v6,v5	#4
		lineVertices[idx] = v2.x;
		lineVertices[idx + 1] = v2.y;
		lineVertices[idx + 2] = v2.z;
		idx += 3;
		lineVertices[idx] = v6.x;
		lineVertices[idx + 1] = v6.y;
		lineVertices[idx + 2] = v6.z;
		idx += 3;
		lineVertices[idx] = v5.x;
		lineVertices[idx + 1] = v5.y;
		lineVertices[idx + 2] = v5.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_up.x;
			trackColors[idxx + 1] = n_up.y;
			trackColors[idxx + 2] = n_up.z;
			idxx += 3;
		}
		// v2,v3,v6	#5
		lineVertices[idx] = v2.x;
		lineVertices[idx + 1] = v2.y;
		lineVertices[idx + 2] = v2.z;
		idx += 3;
		lineVertices[idx] = v3.x;
		lineVertices[idx + 1] = v3.y;
		lineVertices[idx + 2] = v3.z;
		idx += 3;
		lineVertices[idx] = v6.x;
		lineVertices[idx + 1] = v6.y;
		lineVertices[idx + 2] = v6.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_left.x;
			trackColors[idxx + 1] = n_left.y;
			trackColors[idxx + 2] = n_left.z;
			idxx += 3;
		}
		// v3,v7,v6	#6
		lineVertices[idx] = v3.x;
		lineVertices[idx + 1] = v3.y;
		lineVertices[idx + 2] = v3.z;
		idx += 3;
		lineVertices[idx] = v7.x;
		lineVertices[idx + 1] = v7.y;
		lineVertices[idx + 2] = v7.z;
		idx += 3;
		lineVertices[idx] = v6.x;
		lineVertices[idx + 1] = v6.y;
		lineVertices[idx + 2] = v6.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_left.x;
			trackColors[idxx + 1] = n_left.y;
			trackColors[idxx + 2] = n_left.z;
			idxx += 3;
		}
		// v3,v4,v7	#7
		lineVertices[idx] = v3.x;
		lineVertices[idx + 1] = v3.y;
		lineVertices[idx + 2] = v3.z;
		idx += 3;
		lineVertices[idx] = v4.x;
		lineVertices[idx + 1] = v4.y;
		lineVertices[idx + 2] = v4.z;
		idx += 3;
		lineVertices[idx] = v7.x;
		lineVertices[idx + 1] = v7.y;
		lineVertices[idx + 2] = v7.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_down.x;
			trackColors[idxx + 1] = n_down.y;
			trackColors[idxx + 2] = n_down.z;
			idxx += 3;
		}
		// v3,v0,v4	#8
		lineVertices[idx] = v3.x;
		lineVertices[idx + 1] = v3.y;
		lineVertices[idx + 2] = v3.z;
		idx += 3;
		lineVertices[idx] = v0.x;
		lineVertices[idx + 1] = v0.y;
		lineVertices[idx + 2] = v0.z;
		idx += 3;
		lineVertices[idx] = v4.x;
		lineVertices[idx + 1] = v4.y;
		lineVertices[idx + 2] = v4.z;
		idx += 3;
		for (int j = 0; j < 3; ++j) { // color vertices
			trackColors[idxx] = n_down.x;
			trackColors[idxx + 1] = n_down.y;
			trackColors[idxx + 2] = n_down.z;
			idxx += 3;
		}
	}
}

void initGround() {
	// create an integer handle for the texture
	glGenTextures(1, &texHandle);
	int code = initTexture("ground.jpg", texHandle);
	if (code != 0)
	{
		printf("Error loading the texture image.\n");
		exit(EXIT_FAILURE);
	}

	// set up coordinates
	groundVertices = (float*)malloc(sizeof(float) * (3 * 3 * 2));
	groundTexVertices = (float*)malloc(sizeof(float) * (2 * 3 * 2));

	// ground on x-y plane
	glm::vec3 p0 = glm::vec3(-1000.0f, -3000.0f, -100.0f);
	glm::vec3 p1 = glm::vec3(-1000.0f, 3000.0f, -100.0f);
	glm::vec3 p2 = glm::vec3(3000.0f, -1000.0f, -100.0f);
	glm::vec3 p3 = glm::vec3(3000.0f, 3000.0f, -100.0f);

	float amt = 72.0f;
	glm::vec2 v0 = glm::vec2(0.0f, 0.0f);
	glm::vec2 v1 = glm::vec2(0.0f, amt);
	glm::vec2 v2 = glm::vec2(amt, 0.0f);
	glm::vec2 v3 = glm::vec2(amt,amt);

	groundVertices[0] = p0.x; groundVertices[1] = p0.y; groundVertices[2] = p0.z;
	groundVertices[3] = p1.x; groundVertices[4] = p1.y; groundVertices[5] = p1.z;
	groundVertices[6] = p2.x; groundVertices[7] = p2.y; groundVertices[8] = p2.z;

	groundVertices[9] = p0.x; groundVertices[10] = p0.y; groundVertices[11] = p0.z;
	groundVertices[12] = p2.x; groundVertices[13] = p2.y; groundVertices[14] = p2.z;
	groundVertices[15] = p3.x; groundVertices[16] = p3.y; groundVertices[17] = p3.z;

	groundTexVertices[0] = v0.x; groundTexVertices[1] = v0.y;
	groundTexVertices[2] = v1.x; groundTexVertices[3] = v1.y;
	groundTexVertices[4] = v2.x; groundTexVertices[5] = v2.y;

	groundTexVertices[6] = v0.x; groundTexVertices[7] = v0.y;
	groundTexVertices[8] = v2.x; groundTexVertices[9] = v2.y;
	groundTexVertices[10] = v3.x; groundTexVertices[11] = v3.y;

	// creat vbo
	glGenBuffers(1, &groundVertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, groundVertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (3 * 3 * 2 ), groundVertices, GL_STATIC_DRAW);

	glGenBuffers(1, &groundTexVertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, groundTexVertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (2 * 3 * 2), groundTexVertices, GL_STATIC_DRAW);

	textureProgram = new BasicPipelineProgram;
	int ret = textureProgram->Init(textureBasePath);
	if (ret != 0) abort();

	glGenVertexArrays(1, &groundVertexArray);
	glBindVertexArray(groundVertexArray); // bind the VAO

	glBindBuffer(GL_ARRAY_BUFFER, groundVertexBuffer);
	GLuint loc =
		glGetAttribLocation(textureProgram->GetProgramHandle(), "position");
	glEnableVertexAttribArray(loc);
	glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	// bind the VBO “buffer” (must be previously created)
	glBindBuffer(GL_ARRAY_BUFFER, groundTexVertexBuffer);
	// get location index of the “texCoord” shader variable
	loc = glGetAttribLocation(textureProgram->GetProgramHandle(), "texCoord");
	glEnableVertexAttribArray(loc); // enable the “texCoord” attribute
	// set the layout of the “texCoord” attribute data
	glVertexAttribPointer(loc, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindVertexArray(0);
}

void initScene(int argc, char* argv[])
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    // modify the following code accordingly
    int posSize = 0;
    for (int i = 0; i < numSplines; ++i) {
        posSize += splines[i].numControlPoints;
    }
    lineVertices = (float*)malloc(sizeof(float) * 6006 * 24 * posSize);
	//trackColors = (float*)malloc(sizeof(float) * 8008 * 24 * posSize);
	trackColors = (float*)malloc(sizeof(float) * 6006 * 24 * posSize);

	initRailCrossSectionData();

    glGenBuffers(1, &triVertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, triVertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6006 * 24 * posSize, lineVertices,
        GL_STATIC_DRAW);

    glGenBuffers(1, &triColorVertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, triColorVertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6006 * 24 * posSize, trackColors, GL_STATIC_DRAW);

    pipelineProgram = new BasicPipelineProgram;
    int ret = pipelineProgram->Init(shaderBasePath);
    if (ret != 0) abort();

    glGenVertexArrays(1, &triVertexArray);
    glBindVertexArray(triVertexArray);
    glBindBuffer(GL_ARRAY_BUFFER, triVertexBuffer);

    GLuint loc =
        glGetAttribLocation(pipelineProgram->GetProgramHandle(), "position");
    glEnableVertexAttribArray(loc);
    glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, triColorVertexBuffer);
	loc = glGetAttribLocation(pipelineProgram->GetProgramHandle(), "normal");
	glEnableVertexAttribArray(loc);
	glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindVertexArray(0);

	initGround();

    glEnable(GL_DEPTH_TEST);

    sizeTri = 6006 * 24 * posSize;

    std::cout << "GL error: " << glGetError() << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "The arguments are incorrect." << endl;
        cout << "usage: ./hw1 <heightmap file>" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Initializing GLUT..." << endl;
    glutInit(&argc, argv);

    cout << "Initializing OpenGL..." << endl;

    // load the splines from the provided filename
    loadSplines(argv[1]);

    printf("Loaded %d spline(s).\n", numSplines);
    for (int i = 0; i < numSplines; i++)
        printf("Num control points in spline %d: %d.\n", i, splines[i].numControlPoints);


#ifdef __APPLE__
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
#else
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
#endif

    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutCreateWindow(windowTitle);

    cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
    cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
    cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

#ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(windowWidth - 1, windowHeight - 1);
#endif

	initLineData();
	cout << "Loaded line data.\n";
	initCamData();
	cout << "Loaded cam data.\n";

    // tells glut to use a particular display function to redraw 
    glutDisplayFunc(displayFunc);
    // perform animation inside idleFunc
    glutIdleFunc(idleFunc);
    // callback for mouse drags
    glutMotionFunc(mouseMotionDragFunc);
    // callback for idle mouse movement
    glutPassiveMotionFunc(mouseMotionFunc);
    // callback for mouse button changes
    glutMouseFunc(mouseButtonFunc);
    // callback for resizing the window
    glutReshapeFunc(reshapeFunc);
    // callback for pressing the keys on the keyboard
    glutKeyboardFunc(keyboardFunc);

    // init glew
#ifdef __APPLE__
  // nothing is needed on Apple
#else
  // Windows, Linux
    GLint result = glewInit();
    if (result != GLEW_OK)
    {
        cout << "error: " << glewGetErrorString(result) << endl;
        exit(EXIT_FAILURE);
    }
#endif

    // do initialization
    initScene(argc, argv);

    // sink forever into the glut loop
    glutMainLoop();
}


