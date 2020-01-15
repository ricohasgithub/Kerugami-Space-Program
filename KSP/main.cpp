
#include <windows.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

// This struct represents one unit of a point in 3-dimensional space
struct Point3D
{
    double x;
    double y;
    double z;
};

// This struct is used to represent one object (loaded from a .obj file). It contains a list of the points and all the vectors for drawing the faces. It contains the max and min coordinates of every axis for scaling. It also contains physics info about the use of the component.
struct Object
{

    vector<Point3D> vertices;
    vector<Point3D> normals;
    vector<int> triangles;
    vector<int> polygons;
    vector<int> elements;

    double maxX;
    double minX;

    double maxY;
    double minY;

    double maxZ;
    double minZ;

    double mass;
    double thrust;
    double lift;
    double drag;

};

struct Union
{
    vector<vector<Object> > components;
};

// This function returns the magnitude of a given Point3D vector
double getMagnitude (Point3D p)
{
    return sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z));
}

// This method returns a new Point3D with subtracted values
Point3D subtractP3D (Point3D p1, Point3D p2)
{
    Point3D p3;
    p3.x = p1.x - p2.x;
    p3.y = p1.y - p2.y;
    p3.z = p1.z - p2.z;
    return p3;
}

// This method does cross multiplication on the 2 input points
Point3D crossMultiplyP3D (Point3D p1, Point3D p2)
{
    Point3D p3;
    p3.x = (p1.y * p2.z) - (p2.y * p1.z);
    p3.x = (p1.z * p2.x) - (p2.z * p1.x);
    p3.x = (p1.x * p2.y) - (p2.x * p1.y);
    return p3;
}

// This method returns a new normalized unit vector Point3D of given Point3D p
Point3D normalize (Point3D p)
{
    double magnitude = getMagnitude(p);
    Point3D np;
    np.x = (p.x/magnitude);
    np.y = (p.y/magnitude);
    np.z = (p.z/magnitude);
    return np;
}

// This method returns a new loaded .obj file into the program as a new vector of objects
vector<Object> loadObject (string fName) {

    // New vector to load the program files to
    vector<Object> objects;

    // Create a new input file stream reader that will parse the .obj object model file
    ifstream fParser(fName.c_str());

    if (!fParser)
    {
        cout << "Invalid File!" << endl;
        return objects;
    }

    // String variable used to store each line
    string line;

    // Create a new blank object for the first case
    Object newO;
    objects.push_back(newO);

    // Integer variable used to keep track of the current index of the object
    int cObj = 0;

    // These double variables are used to record the maximum x, y and z values of the current object. Used later for scaling and normalization
    double maxX, maxY, maxZ = -1000000;
    // These double variables are used to record the minimum x, y and z values of the current object. Used later for scaling and normalization
    double minX, minY, minZ = 1000000;

    // Start the reading call for the entire file
    while (getline(fParser, line))
    {

        // Check to see if the current line decalres a new object
        if (line.substr(0,2) == "o ")
        {

            // Create a new object and add it to the objects vector. Increment the index of the current object
            Object newO;
            objects.push_back(newO);
            cObj++;

        }
        else if (line.substr(0,2) == "v ")
        {

            // The current line contains a point for a vertex. Create a point and push it to the back of the current object's vertex
            Point3D tempP;

            // Read the current point/vertice
            istringstream sParser(line.substr(2));

            // Parse the x, y and z coordinates into the current object
            sParser >> tempP.x;
            //sParser >> tempP.y;
            //sParser >> tempP.z;

            sParser >> tempP.z;
            sParser >> tempP.y;

            // Check to see if each coordinate is the currently assumed greatest value (maximum); if so, update the current max values
            if (tempP.x > maxX) {
                maxX = tempP.x;
            }

            if (tempP.y > maxY) {
                maxY = tempP.y;
            }

            if (tempP.z > maxZ) {
                maxZ = tempP.z;
            }

            // Check to see if each coordinate is the currently assumed smallest value (minimum); if so, update the current min values
            if (tempP.x < minX) {
                minX = tempP.x;
            }

            if (tempP.y < minY) {
                minY = tempP.y;
            }

            if (tempP.z < minZ) {
                minZ = tempP.z;
            }

            // Add the new point into the current vertices vector in the current object being loaded
            objects[cObj].vertices.push_back(tempP);

        }
        else if (line.substr(0,2) == "f ")
        {

            // Get the amount of vertices and determine whether it is a triangle or a polygon

            // Assume that the current line contains a new face composed of 4 vertex indices only for now
            string as, bs, cs, ds;

            // Tokenize the string line and read each token
            istringstream sParser(line.substr(2));

            // Read each index of the eadges used. Assume that there is a fourth vertice (some file formats do)
            sParser >> as;
            sParser >> bs;
            sParser >> cs;
            sParser >> ds;

            // Read the first 3 vertice indexes first for future drawing refenece
            int a = stoi(as.substr(0,as.find("/")));
            int b = stoi(bs.substr(0,bs.find("/")));
            int c = stoi(cs.substr(0,cs.find("/")));

            // Decrement indices (.obj indices frustratingly start at 1)
            a--;
            b--;
            c--;

            // Check to see if a fourth vertice exists
            if (ds != "") {

                // If the last string token isn't blank, a fourth vertice exists and is added

                int d = stoi(ds.substr(0,ds.find("/")));
                d--;

                objects[cObj].polygons.push_back(a);
                objects[cObj].polygons.push_back(b);
                objects[cObj].polygons.push_back(c);
                objects[cObj].polygons.push_back(d);

            } else {

                // If the last string token is blank, fill triangles instead of polygons
                objects[cObj].triangles.push_back(a);
                objects[cObj].triangles.push_back(b);
                objects[cObj].triangles.push_back(c);

            }

        }

    }

    // Update the max and min coordinate parameters in each object
    for (Object &obj : objects) {

        obj.maxX = maxX;
        obj.minX = minX;

        obj.maxY = maxY;
        obj.minY = minY;

        obj.maxZ = maxZ;
        obj.minZ = minZ;

    }

}

// This method takes in an object and scales it down to an input range
Object scaleObject (Object obj, double nMaxX, double nMinX, double nMaxY, double nMinY, double nMaxZ, double nMinZ) {

        double aspectRatio = (obj.maxY - obj.minY) / (obj.maxX - obj.minX);

        nMaxX *= (1/aspectRatio);
        nMinX *= (1/aspectRatio);

        // Object struct for the new point
        Object nObj = obj;

        vector<Point3D> nVertices;

        for (Point3D p : obj.vertices) {

            // Temporary point variable to be pushed back
            Point3D temp;

            // Scale the x, y and z to be between their respective nMax and nMin values (normalization range)

            temp.x = ( ( (p.x - obj.minX) * (nMaxX - nMinX) ) / (obj.maxX - obj.minX) ) + nMinX;
            temp.y = ( ( (p.y - obj.minY) * (nMaxY - nMinY) ) / (obj.maxY - obj.minY) ) + nMinY;
            temp.z = ( ( (p.z - obj.minZ) * (nMaxZ - nMinZ) ) / (obj.maxZ - obj.minZ) ) + nMinZ;

            nVertices.push_back(temp);

        }

        nObj.vertices = nVertices;

        return nObj;

}

// This void method draws vector of objects with a translation of xpos, ypos, zpos
void drawObject (vector<Object> objects, double xpos, double ypos, double zpos) {

    // Draw every object inside the given objects vector
    for (int m=0; m<objects.size(); m++) {

        // Retrieve the vector of vertices (all points to be drawn)
        vector<Point3D> vertices = objects[m].vertices;

        // Retrieve the list of indexs to draw the triangles and polygonal faces
        vector<int> triangles = objects[m].triangles;
        vector<int> polygons = objects[m].polygons;

        // Draw all triangles
        for (int i=0; i<triangles.size(); i+=3) {

            // First element refers to index of vertex, draw triangle based off of points
            Point3D cVertex1 = vertices[triangles[i]];
            Point3D cVertex2 = vertices[triangles[i+1]];
            Point3D cVertex3 = vertices[triangles[i+2]];

            // Draw each triangle as a new line frame
            glBegin(GL_LINE_STRIP);

                glVertex3d(cVertex1.x + xpos, cVertex1.y + ypos, cVertex1.z + zpos);
                glVertex3d(cVertex2.x + xpos, cVertex2.y + ypos, cVertex2.z + zpos);
                glVertex3d(cVertex3.x + xpos, cVertex3.y + ypos, cVertex3.z + zpos);

            glEnd();

        }

        // Draw all quadrealateral polygonal frames
        for (int i=0; i<polygons.size(); i+=4) {

            // First element refers to index of vertex, draw polygon based off of points
            Point3D cVertex1 = vertices[polygons[i]];
            Point3D cVertex2 = vertices[polygons[i+1]];
            Point3D cVertex3 = vertices[polygons[i+2]];
            Point3D cVertex4 = vertices[polygons[i+3]];

            // Draw the shape as a line frame
            glBegin(GL_LINE_STRIP);

                // Draw the four edges of the line strip (with translation)
                glVertex3d(cVertex1.x + xpos, cVertex1.y + ypos, cVertex1.z + zpos);
                glVertex3d(cVertex2.x + xpos, cVertex2.y + ypos, cVertex2.z + zpos);
                glVertex3d(cVertex3.x + xpos, cVertex3.y + ypos, cVertex3.z + zpos);
                glVertex3d(cVertex4.x + xpos, cVertex4.y + ypos, cVertex4.z + zpos);

            glEnd();

        }

    }
}

// This void method renders a string (s) onto the screen at the given coordinates x, y with a given font
void renderString (double x, double y, void* font, string s) {

    // Set the rasterization coordinates
    glRasterPos2i(x, y);

    // Iterate through each character and render it independently
    for (char c : s) {

        // Draw the bitmap of the current character
        glutBitmapCharacter(font, c);

    }

}

// This integer will represent the current stage of the game
int stage = 0;

// All Rocket objects/components loaded into the program to be used for custom rocket construction
vector<vector<Object> > components;
// The menu displaying the components for the player to select and add to the assembly
vector<vector<Object> > menu;
// This is a list of all objects that a user has selected but not applied to the rocket (i.e., in the "workspace" but not in assembly)
vector<vector<Object> > workspace;

// The union assembly variable used to represent the final assembly to be used in the simulation
Union assembly;

// Index/ID of the currently selected component to be moved (initialize with null value)
int selected = -1;
// Index/ID of the currently selected menu part to be added to the main assembly (initialize with null value)
int menuSelection = -1;

// The x and y coordinates of the left mouse click, scaled between 0 and 1000 (used to determine which menu item has been clicked
double leftX = -1;
double leftY = -1;

// The translation values of the currently selected component
double sxpos = 0;
double sypos = 0;
double szpos = 0;

// Boolean variable stating whether the middle mouse button is being held
bool mhold = false;
// Global Perspective (gp) middle mouse button values (delta x, delta y, current x + y)
int gpx1, gpx2, gpy1, gpy2, gpcx, gpcy;

// Physics engine values - The following variables are constants for simulating a rocket launch

// The current state of the rocket (assembly)
double totalMass = 0;
double totalThrust = 0;
double totalLift = 0;
double totalDrag = 0;

// The vertical position at any point in time
double v_pos = 0.0;
// The vertical velocity at any point in time
double v_vel = 0.0;

// The last time the time is recorded (used to record the change in, delta, time)
double dtime = 0.0;

// Gravitational acceleration on the surface of the Earth (assuming no increase in decelleration as rocket gors further up)
const double g_accl = -9.8;

// This void method takes in a list of items to be displayed "rotating" in display menu. The screen is assumed to be (0, 1000, 0, 1000, -1000, 1000). The method pipes the scaled models to the menu global vector
void initMenu (vector<vector<Object> > items) {

    // Draw each of the objects in the menu. Scale down first and then draw
    for (vector<Object> objList : items) {

        // Create a new temporary vector of objects to store the scaled model
        vector<Object> temp;

        for (Object obj : objList) {

            Object nObj = scaleObject(obj, 200, 0, 200, 0, 200, -200);
            temp.push_back(nObj);

        }

        // Push the new temporary vector into the global menu variable
        menu.push_back(temp);

    }

}

// This void method draws the intro screen
void drawIntroScreen () {

    // Draw a black background
    glClearColor(0.0, 0.0, 0.0, 0.0);

    // Render the introduction text across the screen

    // Set the draw color to white
    glColor3f(1.0, 1.0, 1.0);

    // Render the intro text as seperate lines
    renderString(10, 700, GLUT_BITMAP_HELVETICA_18, "You are the fearless Space Commander in Chief (SCC) of an alien race");
    renderString(10, 650, GLUT_BITMAP_HELVETICA_18, "named the Kerugami. Posed with the previously thought insurmountable");
    renderString(10, 600, GLUT_BITMAP_HELVETICA_18, "objective of making the Kerugami a space-faring species, you are given");
    renderString(10, 550, GLUT_BITMAP_HELVETICA_18, "the necessary components to build and launch a rocket. ");
    renderString(200, 400, GLUT_BITMAP_HELVETICA_18, "Are you ready to take on this daunting task?");
    renderString(670, 250, GLUT_BITMAP_HELVETICA_18, " Press SPACE for YES");
    renderString(520, 200, GLUT_BITMAP_HELVETICA_12, " [No has not been added as a feature to the game]");

}

// This void method draws the menu vector on the side
void drawMenu () {

    // Double variables used to keep track of the starting position of each draw
    double startX = 5.0;
    double startY = 725.0;
    double startZ = 1000.0;

    for (vector<Object> component : menu) {

        // Set the polyon color to white
        glColor3f(1, 1, 1);

        // Draw a square polygon surrounding the current object
        glBegin(GL_POLYGON);
            glVertex3d(startX + 250, startY, startZ);
            glVertex3d(startX + 250, startY + 250, startZ);
            glVertex3d(startX, startY + 250, startZ);
            glVertex3d(startX, startY, startZ);
        glEnd();

        // Draw a bounding box around the polygon

        // Set the line drawing color to black
        glColor3f(0, 0, 0);

        // Draw a square polygon surrounding the current object
        glBegin(GL_LINE_LOOP);
            glVertex3d(startX, startY, 1200);
            glVertex3d(startX + 250, startY, 1200);
            glVertex3d(startX + 250, startY + 250, 1200);
            glVertex3d(startX, startY + 250, 1200);
        glEnd();

        // Draw the actual mini-sized model at the correct starting position

        glColor3f(0, 0, 1);

        drawObject(component, startX, startY, startZ);
        startY -= 250;

    }

    // Draw user instructions (in text) underneath menu

    // Set the text drawing color to black
    glColor3f(0.0, 0.0, 0.0);

    // Draw text with user instruction menu
    renderString(10, 180, GLUT_BITMAP_HELVETICA_12, "Use middle mouse button to rotate");
    renderString(10, 145, GLUT_BITMAP_HELVETICA_12, "Left click on menu item to add a part");
    renderString(10, 110, GLUT_BITMAP_HELVETICA_12, "Press 0-9 to select an unassembled part");
    renderString(10, 75, GLUT_BITMAP_HELVETICA_12, "Press W,A,S,D,P,L to move selected part");
    renderString(10, 40, GLUT_BITMAP_HELVETICA_12, "Press U to assemble wokspace");
    renderString(10, 15, GLUT_BITMAP_HELVETICA_12, "(Assembled parts cannot be moved)");

}

// This method checks to see if any components needs to be added to the workspace
void updateWorkspace () {// Scale the creen for the specific screen

    // Check to see if the mouse left click selection lands on a valid menu item. Each menu item is bounded by a 250 by 250 box
    if (leftX <= 250) {

        // Get the current index of the selected component in the workspace based off the click position
        int index = (int) (leftY/250);

        // Check to see if the X value is inside the range of the menu
        if (index < menu.size()) {
            // Update the menu item at the selected menu "square" by adding it to the workspace
            workspace.push_back(components[index]);
        }

    }

}

// This void method takes in a nested vector of objects and pipes all subobjects as objects in the assembly union. It also clears the entire workspace
void assembleComponents(vector <vector <Object> > parts) {

    // Iterate through every object group
    for (int i=0; i<parts.size(); i++) {

        // Get the current set of components to be added to the assembly and removed from the workspace
        vector<Object> objects = parts[i];

        // Add the current sub-component/object group to the main assembly
        assembly.components.push_back(objects);
        workspace.erase(workspace.begin() + i);

    }

}

// This void method increments all the points at index in the components vector (pre-setting a translation)
void setPreTranslate (vector< vector <Object> > &components, int index, int nx, int ny, int nz) {

    // Get the pointer objects vector to be modified
    vector<Object> &objects = components[selected];

    // Iterate through all "sub-objects" within the current object and update each vertice seperately
    for (int m=0; m<objects.size(); m++) {

        // Get all vertices in the current sub-object
        vector<Point3D> &vertices = objects[m].vertices;

        // Iterate over all of the current vertices and increment each
        for (int v=0; v<vertices.size(); v++) {
            // Increment the current vertice
            Point3D &cVertex = vertices[v];
            cVertex.x += nx;
            cVertex.y += ny;
            cVertex.z += nz;
        }

    }

}

// This void method draws the entire rocket assembly screen
void drawRocketAssembly () {

    // Draw a white background
    glClearColor(1.0, 1.0, 1.0, 0.0);

    // Draw each of the different components in the current assembly
    for (int i=0; i<assembly.components.size(); i++) {

        // Get the current obect to be drawn
        vector<Object> obj = assembly.components[i];

        double xtrans = 500;
        double ytrans = 500;
        double ztrans = 0;

        // Set color to the completed assembly union color black
        glColor3d(0,0,0);

        glPushMatrix();

            glMatrixMode(GL_MODELVIEW);

            // Rotate the global perspective
            glRotated(gpcx, 0, 1000, 0);
            glRotated(gpcy, 1000, 0, 0);

            drawObject(obj, xtrans, ytrans, ztrans);

        glPopMatrix();

    }

    // Draw each of the different components in the current workspace (excluding assembly)
    for (int i=0; i<workspace.size(); i++) {

        // Get the current obect to be drawn
        vector<Object> obj = workspace[i];

        double xtrans = 500;
        double ytrans = 500;
        double ztrans = 0;

        // Determine drawing color depending on whether the current element is the selected one to be moved
        if (i == selected) {
            // The element is selected, set draw color to green
            glColor3d(0,1,0);
            // Also move the component
            xtrans += sxpos;
            ytrans += sypos;
            ztrans += szpos;

        } else {
            // Otherwise, not selected. Set to default draw color (blue)
            glColor3d(0,0,1);
        }

        glPushMatrix();

            glMatrixMode(GL_MODELVIEW);

            // Rotate the global perspective
            glRotated(gpcx, 0, 1000, 0);
            glRotated(gpcy, 1000, 0, 0);

            drawObject(obj, xtrans, ytrans, ztrans);

        glPopMatrix();

    }

    // Draw the components menu
    drawMenu();

}

// This void method sets up the constants for the current iteration of the rocket (assembly)
void updateRocketPhysics () {

    // Zero all previous rocket values
    totalMass = 0;
    totalThrust = 0;
    totalLift = 0;
    totalDrag = 0;

    // Get the total values for mass, thrust, lift and drag on the rocket (assembly). Iterate through the entire assembly and retireve all components (and included subcomponent) data
    for (vector<Object> temp_components : assembly.components) {

        // Iterate through all sub components in the assembly
        for (Object obj : temp_components) {

            // Accumulate physics engine values
            totalMass += obj.mass;
            totalThrust += obj.thrust;
            totalLift += obj.lift;
            totalDrag += obj.drag;

        }

    }

    // Update the initial vertical velocity to the total starting thrust
    v_vel = totalThrust;

}

// This void method does translation operations on the assembled rocket based on its physics engine values
void launchRocket () {

    if (v_pos < 0) {

        // If the rocket crashed (goes below groun), stop further calls
        return;

    }

    // Get the curretn elapsed time
    const double time = glutGet(GLUT_ELAPSED_TIME) / 10000.0;

    // The current vertical acceleration
    double v_accl = g_accl;

    // Calculate the horizontal and vertical acceleration and velocity: Thrust is the initial vertical velocity, Lift is the initial vertical accelleration, and Drag is the rate of vertical decelleration (for thrust to decrease till 0)

    // First, calculate additional vertical acceleration (lift minus drag)
    if (totalLift > 0) {

        // Get the total vertical acceleration by adding the current total lift to it
        v_accl += totalLift;

        // Decrease lift further if it is not already zero
        totalLift -= totalDrag;

    } else if (totalLift < 0) {

        // If the lift is less than 0, set to to zero (drag should only decrease additional vertical acceleration)
        totalLift = 0;

        // Reset the vertical acceleration to only the gravitational decelleration
        v_accl = g_accl;

    }

    // Get the change in time
    dtime = time - dtime;

    // Update the vertical velocity (accumulated acceleration over time)
    v_vel += v_accl * dtime;
    // update the vertical position (accumulated velocity over time)
    v_pos += v_vel * dtime;

    cout << "Position: " << v_pos << endl;
    cout << "Velocity: " << v_vel << endl;
    cout << "Acceleration: " << v_accl << endl;



    // Draw the rocket

}

// This void method draws the rocket launch simulation
void drawRocketLaunch () {

    // Draw a white background
    glClearColor(1.0, 1.0, 1.0, 0.0);

    // Draw each of the different components in the current assembly
    for (int i=0; i<assembly.components.size(); i++) {

        // Get the current obect to be drawn
        vector<Object> obj = assembly.components[i];

        // Set color to the completed assembly union color black
        glColor3d(0,0,0);

        glPushMatrix();

            // Initialize the matrix modelview mode for glDraw display
            glMatrixMode(GL_MODELVIEW);

            // Draw the rocket
            drawObject(obj, 500, v_pos, 0);

        glPopMatrix();

    }
}

// This is the default display method called by redraws. It contains code to distinguish the current stage of the game and draw the appropriate screen
void display(void) {

    // Clear the current color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Draw the appropriate "stage" of the game depending on the user response
    switch (stage) {

        // Stage 0: intro screen
        case 0:

            // Scale the creen for the specific screen
            glLoadIdentity();
            glOrtho(0, 1000, 0, 1000, -1200, 1200);

            // Draw the intro screen
            drawIntroScreen();
            break;

        // Stage 1: rocket assembly screen
        case 1:

            // Scale the creen for the specific screen
            glLoadIdentity();
            glOrtho(0, 1000, 0, 1000, -1200, 1200);

            // Draw the rocket assembly screen
            drawRocketAssembly();
            break;

        // Stage 2: rocket launch screen
        case 2:

            // Scale the creen for the specific screen
            glLoadIdentity();
            glOrtho(0, 20000, 0, 20000, -1200, 1200);

            // Run rocket launching simulation
            launchRocket();
            // Draw the rocket once simulation phase completes for current cycle
            drawRocketLaunch();

            break;

    }

    // Flushing the matrix to the cache for double buffering
    glutSwapBuffers();

}

// Set the idle animation
void idle(void) {
    glutPostRedisplay();
}

// This void method takes in a filepath/filename for the components text file and then buffers and prepares the entire components vector
void loadComponents (string filename) {

    // Initialize a new file parser to read from the components filename
    ifstream fParser(filename.c_str());

    if (!fParser)
    {
        cout << "Invalid File!" << endl;
        return;
    }

    // Temporary string variable used to read every component filename at every line
    string componentFileName;

    // Read through every line in the file and load every file into the components vector
    while (getline(fParser, componentFileName)) {

        // Load the current line/filename into the components vector as a new object
        vector<Object> temp_components = loadObject(componentFileName);

        // Temporary string variable used to read the following lines for physics engine data
        string temp;

        // Read the physics data from the following lines and add them to the component. The order goes: mass, thrust, lift and drag
        getline(fParser, temp);
        double mass = stod(temp);
        getline(fParser, temp);
        double thrust = stod(temp);
        getline(fParser, temp);
        double lift = stod(temp);
        getline(fParser, temp);
        double drag = stod(temp);

        // Add the data to each sub-component of the component
        for (int i=0; i<temp_components.size(); i++) {

            // Get the current component update its physics engine parameters
            temp_components[i].mass = mass;
            temp_components[i].thrust = thrust;
            temp_components[i].lift = lift;
            temp_components[i].drag = drag;

        }

        // Add the component into the components vector
        components.push_back(temp_components);

    }

}

// This method opens the .obj files and parses them to load the object models
void init() {

    // The number of components
    int fLength = 3;

    // Load all the components into the components vector
    loadComponents("C://Users/ricoz/Desktop/C++ Workspace/Kerugami-Space-Program/KSP/Components.txt");

    // Initialize the menu once the components have been loaded
    initMenu(components);

}

// This method enables global perspective via the middle mouse button (during rocket assembly stage)
void manipulateObjects (int x, int y) {

    if (stage == 1 && mhold) {

        // Center mouse button is being held while in the rocket assembly stage, update newest x and y coordinate changes

        gpx1 = gpx2;
        gpy1 = gpy2;

        gpx2 = x;
        gpy2 = y;

        // Get the change in x and y
        double dx = gpx2 - gpx1;
        double dy = gpy2 - gpy1;

        // Update cx and cy (the current rotation angle) depending on whether dx and dy was increasing or decreasing
        if (dx < 0) {
            // Move to the left
            gpcx += 1;
        } else if (dx > 0) {
            // Move the the right
            gpcx -= 1;
        }

        if (dy < 0) {
            // Move up
            gpcy -= 1;
        } else if (dy > 0) {
            // Move down
            gpcy += 1;
        }

        // Update the screen
        glutPostRedisplay();

    }

}

// This method will listen for all mouse button controls. This includes adjusting the global assembly perspective and
void mouseListner (int button, int state, int x, int y) {

    // Global viewing perspective changes (middle mouse button) in the rocket assembly stage (stage = 1)

    // Get the change in x and change in x to determine rotation angle. Use solidworks controls (center mouse button while held)
    if (stage == 1 && button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {

        // Update mouse down to alert mouse move check
        mhold = true;

        // Center mouse button is being held, update newest x and y coordinate changes

        gpx1 = gpx2;
        gpy1 = gpy2;

        gpx2 = x;
        gpy2 = y;

    } else {
        // Mouse button is not being held
        mhold = false;
    }

    // Menu component selection (left mouse button click) in the rocket assembly stage (stage = 1)

    // Listner activated for a menu selection. Update the click position
    if (stage == 1 && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {

        // Scale the click position to be between 0 and 1000 instead of 0 and 600

        leftX = (double) x * (1000.0/600.0);
        leftY = (double) y * (1000.0/600.0);

        // Update the workspace
        updateWorkspace();

    }

}

// This method listens for keyboard events during the individual stages of the game
void keyboardListener (unsigned char key, int x, int y) {

    // Check for the current stage of the game; key actions will change depending on the state
    switch (stage) {

        case 0:

            // Stage 0: Intro text screen

            if (key == 32) {

                // Player pressed space bar, move onto next stage (rocket assembly, stage 1)
                stage = 1;
                // Draw the new screen
                glutPostRedisplay();

                // Stop the keyboard listening call right away
                break;

            }

        case 1:

            // Stage 1: Rocket assembly stage

            if (key == 32) {

                // Player pressed space bar, move onto next stage (rocket launch, stage 2)
                stage = 2;

                // Update the current state of the newly assembled rocket
                updateRocketPhysics();

                // Draw the new screen
                glutPostRedisplay();

                // Stop the keyboard listening call right away
                break;

            }

            // Check to see if the key pressed was to select a certain component in the workspace (ascii values 48 - 57). Only actiavted while in the rocket assembly stage (stage = 1)
            if (key >= 48 && key <= 57 && (key - 48) < workspace.size()) {

                if (selected != -1) {
                    // Update all the point values of the previous object
                    setPreTranslate(workspace, selected, sxpos, sypos, szpos);
                }

                // The key is a number key and the value is valid (there are that many componenets)
                selected = key - 48;

                // Reset the translation position of the object
                sxpos = 0;
                sypos = 0;
                szpos = 0;

            }

            // Move the appropriate objects if an object is selected and a key is pressed to mvoe the selected component
            if (selected != -1 && key == 'w') {
                // Increase the y value of the selected component
                sypos += 5;
            } else if (selected != -1 && key == 's') {
                // Decrease the y value of the selected component
                sypos -= 5;
            } else if (selected != -1 && key == 'a') {
                // Decrease the x value of the selected component
                sxpos -= 5;
            } else if (selected != -1 && key == 'd') {
                // Increase the x value of the selected component
                sxpos += 5;
            } else if (selected != -1 && key == 'p') {
                // Increase the z value of the selected component
                szpos += 5;
            } else if (selected != -1 && key == 'l') {
                // Decrease the z value of the selected component
                szpos -= 5;
            }

            if (key == 'u') {

                // If there has been a translation/modification
                if (selected != -1) {
                    // Update all the point values of the current workspace
                    setPreTranslate(workspace, selected, sxpos, sypos, szpos);
                }

                // Assemble all existing components in the workspace together
                assembleComponents(workspace);
                // Reset the selected variable as no item is selected
                selected = -1;
            }

            break;

    }

}

// Arrays mapping to texture and lighting/shading constants
const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 1000.0f, 1000.0f, 0.0f, 0.0f };

const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };

int main( int argc, char **argv )
{
    // Initialize the new frame and clear the depth buffer
    glutInit( &argc, argv );
    init();
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize( 600, 600 );
    glutCreateWindow( "GLUT .obj Demo" );

    // Set the display function to draw the solid
    glutDisplayFunc(display);
    // Set the idle animation funciton
    glutIdleFunc(idle);
    // Set the mouse event animation funciton
    glutMouseFunc(mouseListner);
    // Set the mouse motion/move function
    glutMotionFunc(manipulateObjects);
    // Set the keyboard function
    glutKeyboardFunc(keyboardListener);
    // Clear the background
    glClearColor(1,1,1,1);

    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Enable materials rendering
 /*   glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);*/

    glutMainLoop();
    return 0;
}
