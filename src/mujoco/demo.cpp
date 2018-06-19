#include <mujoco.h>
#include <glfw3.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <sstream>

// model
mjModel* m = 0;
mjData* d = 0;

// abstract visualization
mjvCamera cam;                      // abstract camera
mjvPerturb pert;                    // perturbation object
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

using namespace std; 

// simple controller applying damping to each dof
void mycontroller(const mjModel* m, mjData* d)
{
    if( m->nu==m->nv )
        mju_scl(d->ctrl, d->qvel, -0.1, m->nv);
}

// get vector of times that correspond to controls
bool read_times(string filepath, vector<double> &times) {

   // open file
   ifstream in(filepath);

   // Check if object is valid
   if(!in)
   {
      cerr << "Cannot open the file : "<< filepath << endl;
      return false;
   }
 
   string str;
   // Read the next line from File untill it reaches the end.
   while (in.good())
   {
      string substr;
      getline(in, substr, ',');
      double t = atof(substr.c_str());
      times.push_back(t);
   }
   //Close The File
   in.close();
   return true;
}

// get vector of controls
// TODO THIS IS HARDCODED FOR U_DIM = 2
bool read_controls(string filepath, vector<vector<double> > &controls){
   // open file
   ifstream in(filepath);

   // Check if object is valid
   if(!in)
   {
      cerr << "Cannot open the file : "<< filepath <<endl;
      return false;
   }
 
   string str;
   // line 1 is u1(t), line 2 is u2(t), ...
   int ctrl_idx = 0;
   string line;
   while (std::getline(in, line))
   {
      stringstream ss(line);
      string substr;
      int t = 0;
      while(getline(ss, substr, ',')){
         double u = atof(substr.c_str());
         // first iteration through, need to create each interior vector
         if(ctrl_idx == 0){
            vector<double> uvec;
            uvec.push_back(u);
            controls.push_back(uvec);
         }else{
            controls[t].push_back(u);
         }
         t++;
      }
      ctrl_idx++;
   }
   //Close The File
   in.close();
   return true;
}

// get vector of states
// TODO THIS IS HARDCODED FOR U_DIM = 2
bool read_states(string filepath, vector<vector<double> > &states){
   // open file
   ifstream in(filepath);

   // Check if object is valid
   if(!in)
   {
      cerr << "Cannot open the file : "<< filepath <<endl;
      return false;
   }
 
   string str;
   // line 1 is u1(t), line 2 is u2(t), ...
   int ctrl_idx = 0;
   string line;
   while (std::getline(in, line))
   {
      stringstream ss(line);
      string substr;
      int t = 0;
      while(getline(ss, substr, ',')){
         double x = atof(substr.c_str());
         // first iteration through, need to create each interior vector
         if(ctrl_idx == 0){
            vector<double> xvec;
            xvec.push_back(x);
            states.push_back(xvec);
         }else{
            states[t].push_back(x);
         }
         t++;
      }
      ctrl_idx++;
   }
   //Close The File
   in.close();
   return true;
}

// for debugging
void print_times(vector<double> times){
   cout << "-------------\n";
   cout << "times: ";
   for(int idx = 0; idx<times.size(); idx++){
      cout << times[idx] << "\n";
   }
   cout << "-------------\n";
}

void print_controls(vector<vector<double> > controls){
   cout << "-------------\n";
   cout << "controls [u1, u2]: \n";
   for(int idx = 0; idx<controls.size(); idx++){
      cout << controls[idx][0] << ", " << controls[idx][1] << "\n";
   }
   cout << "-------------\n";
}

void print_states(vector<vector<double> > states){
   cout << "-------------\n";
   cout << "states [x, theta]: \n";
   for(int idx = 0; idx<states.size(); idx++){
      cout << states[idx][0] << ", " << states[idx][1] << "\n";
   }
   cout << "-------------\n";
}

// interpolates control signal based on current time
void interpolate_controls(double t, vector<double> times, vector<vector<double> > controls, vector<double> &u){
   // sanity check
   u.clear();

   int tf_idx = times.size()-1;
   double tf = times[tf_idx];
   cout << "curr time: " << t << ", final time: " << tf << endl;
   if(t > tf){
      u = controls[tf_idx];
   }else{
      int curr_idx = lower_bound(times.begin(), times.end(), t) - times.begin();
      if(curr_idx == tf_idx){
          u = controls[tf_idx];
          return;
      }
      vector<double> prev = controls[curr_idx];
      vector<double> next = controls[curr_idx+1];
      double ti = times[curr_idx];
      double tf = times[curr_idx+1];
      double u1 = (next[0] - prev[0])*((t-ti)/(tf-ti)) + prev[0];
      double u2 = (next[1] - prev[1])*((t-ti)/(tf-ti)) + prev[1];
      u.push_back(u1);
      u.push_back(u2);
   }
}

// interpolates state based on current time
void interpolate_states(double t, vector<double> times, vector<vector<double> >  states, vector<double> &x){
   // sanity check
   x.clear();

   int tf_idx = times.size()-1;
   double tf = times[tf_idx];
   cout << "curr time: " << t << ", final time: " << tf << endl;
   if(t > tf){
      x = states[tf_idx];
   }else{
      int curr_idx = lower_bound(times.begin(), times.end(), t) - times.begin();
      if(curr_idx == tf_idx){
          x = states[tf_idx];
          return;
      }
      vector<double> prev = states[curr_idx];
      vector<double> next = states[curr_idx+1];
      double ti = times[curr_idx];
      double tf = times[curr_idx+1];
      double x1 = (next[0] - prev[0])*((t-ti)/(tf-ti)) + prev[0];
      double x2 = (next[1] - prev[1])*((t-ti)/(tf-ti)) + prev[1];
      x.push_back(x1);
      x.push_back(x2);
   }
}


// applies the control from vector u to the mujoco sim
void apply_controls(mjData* d, vector<double> &u){
   d->ctrl[0] = u[0];
   d->ctrl[1] = u[1];
}

int main(void)
{

   string MUJOCO_PATH = "/home/abajcsy/Documents/mujoco/mjpro150";
   string license = MUJOCO_PATH + "/bin/mjkey.txt";
   string times_path = "times.csv";
   string controls_path = "controls.csv";
   string states_path = "states.csv";

   // get controller and times
   vector<double> times;
   vector<vector<double> > controls;
   vector<vector<double> > states;
   read_times(times_path, times);
   read_controls(controls_path, controls);
   read_states(states_path, states);

   print_times(times);
   print_controls(controls);
   print_states(states);

   // activate MuJoCo Pro
   mj_activate(license.c_str());

   char error[1000];
   // load model from file and check for errors: cartpole1link_mujoco.xml
   m = mj_loadXML("../../resources/cartpole1link_mujoco.xml", NULL, error, 1000);
   if( !m )
   {
      printf("%s\n", error);
      return 1;
   }

   // make data corresponding to model
   d = mj_makeData(m);

   glfwInit();
   GLFWwindow* window = glfwCreateWindow(1800, 900, "Demo", NULL, NULL);
   // make context current, disable v-sync
   glfwMakeContextCurrent(window);
   glfwSwapInterval(1);

   // initialize visualization data structures
   mjv_defaultCamera(&cam);
   mjv_defaultPerturb(&pert);
   mjv_defaultOption(&opt);
   mjr_defaultContext(&con);
   mjv_makeScene(&scn, 1000);                     // space for 1000 objects
   mjr_makeContext(m, &con, mjFONTSCALE_100);     // model-specific context

   // get framebuffer viewport
   mjrRect viewport = {0, 0, 0, 0};
   glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

   // install control callback
   //mjfGeneric mjcb_control = mycontroller;

   std::vector<double> u, x;
   mjtNum start_t = d->time;
   cout << "start time: " << start_t << endl; 

   // install control callback
   mjcb_control = mycontroller;

   // run main loop, target real-time simulation and 60 fps rendering
   while( !glfwWindowShouldClose(window))
   {
      mjtNum simstart = d->time;
      // advance interactive simulation for 1/60 sec
      //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
      //  this loop will finish on time for the next frame to be rendered at 60 fps.
      //  Otherwise add a cpu timer and exit this loop when it is time to render.
      while( d->time - simstart < 1.0/60.0 ){
         mj_step1(m, d);

         cout << "d->time: " << d->time << endl; 
         mjtNum curr_t = d->time - start_t;

         //interpolate_controls(curr_t, times, controls, u);
         //apply_controls(d, u);
         interpolate_states(curr_t, times, states, x);
         d->qpos[0] = x[0];
         d->qpos[1] = x[1];

         cout << "control: " << d->ctrl[0] << ", " << d->ctrl[1] << endl;
         mj_step2(m, d);
      }

      // get framebuffer viewport
      mjrRect viewport = {0, 0, 0, 0};
      glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

      // update scene and render
      mjv_updateScene(m, d, &opt, &pert, &cam, mjCAT_ALL, &scn);
      mjr_render(viewport, &scn, &con);

      // swap OpenGL buffers (blocking call due to v-sync)
      glfwSwapBuffers(window);

      // process pending GUI events, call GLFW callbacks
      glfwPollEvents();
   }

   // close GLFW, free visualization storage
   glfwTerminate();
   mjv_freeScene(&scn);
   mjr_freeContext(&con);

   // free model and data, deactivate
   mj_deleteData(d);
   mj_deleteModel(m);
   mj_deactivate();

   return 0;
}
