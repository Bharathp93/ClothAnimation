#include<GL/freeglut.h>
#include<vector>
#include<glm/glm.hpp>
#include<iostream>

using namespace std;

class Spring {
public:
	int p1, p2;				//two points of the spring
	float rest_length;		//rest length
	float Ks;				//spring stiffness
	float Kd;				//damping coefficient
	Spring(int p1, int p2, float rest_length, float Ks, float Kd) : p1{ p1 }, p2{ p2 }, rest_length{ rest_length }, Ks{ Ks }, Kd{ Kd } {};
};

vector<Spring> springs;

vector<glm::vec3>index;
vector<glm::vec3>points;
vector<glm::vec3>velocity;
vector<glm::vec3>force;

int cols = 50, rows = 50;
const int numOfPointMasses = (cols + 1)*(rows + 1);
int csize = 4;
float hsize = csize / 2.0f;
float rot_x = 5, rot_y = 0;
float pos_z = -15;

float timeStep = 1 / 60.0f;
double counter = timeStep;

const float DEFAULT_DAMPING = -0.015f;
float	KsStruct = 0.75f, KdStruct = -0.25f;
float	KsShear = 0.75f, KdShear = -0.25f;
float	KsBend = 0.95f, KdBend = -0.25f;
glm::vec3 gravity = glm::vec3(0.0f, -0.00981f, 0.0f);

float mass = 0.5f;

GLint viewport[4];
GLdouble _matrix[16];
GLdouble projection_matrix[16];

/* Mouse Interface  */
int _mouseX = 0;
int _mouseY = 0;
bool _mouseLeft = false;

int spring_count = 0;

void AddSpring(int a, int b, float ks, float kd) {
	Spring spring = Spring(a, b, sqrt(glm::dot((points[a] - points[b]), (points[a] - points[b]))), ks, kd);
	springs.push_back(spring);
}

void mouseEvent(int button, int state, int x, int y)
{
	_mouseX = x;
	_mouseY = y;

	if (state == GLUT_UP)
		switch (button) {
		case GLUT_LEFT_BUTTON:
			_mouseLeft = false;
			break;
		}
	else
		switch (button) {
		case GLUT_LEFT_BUTTON:
			_mouseLeft = true;
			break;
		default:
			break;
		}
}

void mouseMoveEvent(int x, int y)
{
	const int dx = x - _mouseX;
	const int dy = y - _mouseY;

	if (dx == 0 && dy == 0)
		return;

	else if (_mouseLeft) {

		rot_x += dy;
		rot_y += dx;
	}

	_mouseX = x;
	_mouseY = y;

	glutPostRedisplay();
}

void initialise() {

	glEnable(GL_DEPTH_TEST);
	int v = rows + 1;
	int u = cols + 1;

	force.resize(numOfPointMasses);

	//fill in points and velocity
	for (int j = 0; j < v; j++) {
		for (int i = 0; i < u; i++) {
			points.push_back(glm::vec3(((float(i) / (u - 1)) * 2 - 1)* hsize, csize + 1, ((float(j) / (v - 1))* csize)));
			velocity.push_back(glm::vec3(0, 0, 0));
		}
	}

	//fill in index values
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int a = i * (cols + 1) + j;
			int b = a + 1;
			int c = a + (cols + 1);
			int d = c + 1;
			if ((j + i) % 2) {
				index.push_back(glm::vec3(a, c, b));
				index.push_back(glm::vec3(b, c, d));
			}
			else {
				index.push_back(glm::vec3(a, c, d));
				index.push_back(glm::vec3(a, d, b));
			}
		}
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(2);

	//setup springs
	for (int x = 0; x <= rows; x++)
		for (int y = 0; y < cols; y++) {
			AddSpring((u * x) + y, ((u * x) + y) + 1, KsStruct, KdStruct);
		}

	for (int x = 0; x < rows; x++)
		for (int y = 0; y <= cols; y++) {
			AddSpring((u * x) + y, (u * (x + 1)) + y, KsStruct, KdStruct);
		}

	for (int x = 0; x < rows; x++)
		for (int y = 0; y < cols; y++) {
			AddSpring((u * x) + y, (u * (x + 1)) + (y + 1), KsShear, KdShear);
			AddSpring((u * x) + (y + 1), (u * (x + 1)) + y, KsShear, KdShear);
		}

	for (int x = 0; x <= rows; x++)
		for (int y = 0; y <= cols; y++) {
			if ((x < (rows - 1)) && (y < (cols - 1))) {
				AddSpring((u * x) + y, (u * x) + (y + 2), KsBend, KdBend);
				AddSpring((u * x) + y, (u * (x + 2)) + y, KsBend, KdBend);
				AddSpring((u * x) + y, (u * (x + 2)) + (y + 2), KsBend, KdBend);
				AddSpring((u * x) + (y + 2), (u * (x + 2)) + y, KsBend, KdBend);
			}
			else if ((x < (rows - 1)) && !(y < (cols - 1))) {
				AddSpring((u * x) + y, (u * (x + 2)) + y, KsBend, KdBend);
			}
			else if (!(x < (rows - 1)) && (y < (cols - 1))) {
				AddSpring((u * x) + y, (u * x) + (y + 2), KsBend, KdBend);
			}
		}
}

void OnReshape(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)width / (GLfloat)height, 1.f, 100.0f);
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glMatrixMode(GL_MODELVIEW);
}

void OnRender() {

	const int grid = 10;
	counter += (float)1 / 6;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslatef(0, 0, pos_z);
	glRotatef(rot_x, 1, 0, 0);
	glRotatef(rot_y, 0, 1, 0);

	glBegin(GL_LINES);
	glColor3f(0.5f, 0.5f, 0.5f);
	for (int i = -grid; i <= grid; i++)
	{
		glVertex3f(i, 0, -grid);
		glVertex3f(i, 0, grid);
		glVertex3f(-grid, 0, i);
		glVertex3f(grid, 0, i);
	}
	glEnd();

	glColor3f(0.5, 0.5, 0.5);
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < index.size(); i++) {
		glm::vec3 p1 = points[index[i].x];
		glm::vec3 p2 = points[index[i].y];
		glm::vec3 p3 = points[index[i].z];
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();

	glBegin(GL_POINTS);
	for (size_t i = 0; i < numOfPointMasses; i++) {
		glm::vec3 p = points[i];
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();

	glutSwapBuffers();
}

void calculateForcesOnParticles() {
	size_t i = 0;
	for (i = 0; i < numOfPointMasses; i++) {
		force[i] = glm::vec3(0, 0, 0);

		if (i != 0 && i != cols)
			force[i] += gravity;

		force[i] += DEFAULT_DAMPING * velocity[i];
	}

	//add spring forces
	for (i = 0; i < springs.size(); i++) {
		glm::vec3 deltaP = points[springs[i].p1] - points[springs[i].p2];
		glm::vec3 deltaV = velocity[springs[i].p1] - velocity[springs[i].p2];
		float dist = glm::length(deltaP);

		glm::vec3 springForce = ((-springs[i].Ks * (dist - springs[i].rest_length)) + (springs[i].Kd * (glm::dot(deltaV, deltaP) / dist))) * (deltaP / dist);

		if (springs[i].p1 != 0 && springs[i].p1 != cols)
			force[springs[i].p1] += springForce;
		if (springs[i].p2 != 0 && springs[i].p2 != cols)
			force[springs[i].p2] -= springForce;
	}
}

void IntegrateEuler(float dt) {
	float deltaTimeByMass = dt / mass;

	for (int i = 0; i < numOfPointMasses; i++) {
		glm::vec3 oldV = velocity[i];
		velocity[i] += (force[i] * deltaTimeByMass);
		points[i] += dt * oldV;

		//So that the cloth stays above the platform
		if (points[i].y < 0) {
			points[i].y = 0;
		}
	}

}

void DynamicInverseConstraint() {

	for (size_t i = 0; i < springs.size(); i++) {
		//check the current lengths of all springs
		glm::vec3 deltaP = points[springs[i].p1] - points[springs[i].p2];
		float dist = glm::length(deltaP);
		if (dist > springs[i].rest_length) {
			dist -= (springs[i].rest_length);
			dist /= 2.0f;
			deltaP = glm::normalize(deltaP);
			deltaP *= dist;
			if (springs[i].p1 == 0 || springs[i].p1 == cols) {
				velocity[springs[i].p2] += deltaP;
			}
			else if (springs[i].p2 == 0 || springs[i].p2 == cols) {
				velocity[springs[i].p1] -= deltaP;
			}
			else {
				velocity[springs[i].p1] -= deltaP;
				velocity[springs[i].p2] += deltaP;
			}
		}
	}
}

void OnIdle() {

	if (counter >= timeStep)
	{
		calculateForcesOnParticles();
		IntegrateEuler(timeStep);
		DynamicInverseConstraint();
		counter -= timeStep;
	}
	glutPostRedisplay();
}

void main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(1024, 1024);
	glutCreateWindow("Cloth Animation");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(mouseEvent);
	glutMotionFunc(mouseMoveEvent);

	initialise();

	glutMainLoop();
}
