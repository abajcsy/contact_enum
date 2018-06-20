#include <mujoco.h>
#include <glfw3.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <cstring>

// model
mjModel* m = 0;
mjData* d = 0;

// abstract visualization
mjvCamera cam;                      // abstract camera
mjvPerturb pert;                    // perturbation object
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context
mjvFigure figconstraint;
mjvFigure figcost;
mjvFigure figtimer;
mjvFigure figsize;
mjvOption vopt;
char status[1000] = "";

// OpenGL rendering
int refreshrate;
float depth_buffer[5120*2880];        // big enough for 5K screen
unsigned char depth_rgb[1280*720*3];  // 1/4th of screen

// user state
bool paused = false;
bool showoption = false;
bool showinfo = false;
bool slowmotion = false;
bool showdepth = false;
int showhelp = 0;                   // 0: none; 1: brief; 2: full
bool showprofiler = false;
int keyreset = -1;                  // non-negative: reset to keyframe

// help strings
const char help_title[] = 
"Help\n"
"Option\n"
"Info\n"
"Depth\n"
"Full screen\n"
"Stereo\n"
"Profiler\n"
"Slow motion\n"
"Key reset\n"
"Pause\n"
"Reset\n"
"Forward\n"
"Back\n"
"Forward 100\n"
"Back 100\n"
"Autoscale\n"
"Reload\n"
"Geoms\n"
"Sites\n"
"Select\n"
"Center\n"
"Track\n"
"Zoom\n"
"Translate\n"
"Rotate\n"
"Perturb\n"
"Free Camera\n"
"Camera\n"
"Frame\n"
"Label\n"
"Fontsize";


const char help_content[] = 
"F1\n"
"F2\n"
"F3\n"
"F4\n"
"F5\n"
"F6\n"
"F7\n"
"Enter\n"
"Page Up/Down\n"
"Space\n"
"BackSpace\n"
"Right arrow\n"
"Left arrow\n"
"Down arrow\n"
"Up arrow\n"
"Ctrl A\n"
"Ctrl L\n"
"0 - 4\n"
"Shift 0 - 4\n"
"L dblclick\n"
"R dblclick\n"
"Ctrl R dblclick\n"
"Scroll or M drag\n"
"[Shift] R drag\n"
"L drag\n"
"Ctrl [Shift] L/R drag\n"
"Esc\n"
"[ ]\n"
"; '\n"
". /\n"
"- =";

char opt_title[1000] = "";
char opt_content[1000];

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

// writes a vector of states to file
void write_states(string filename, vector<vector<double> > states){
   ofstream myfile;
   myfile.open (filename);
   int num_states = states[0].size();
   for(int s=0;s<num_states;s++){
      for(int t=0;t<states.size();t++){
         if (t == states.size()-1){
            myfile << states[t][s];  
         }else{
            myfile << states[t][s] << ", ";
         }
      } 
      myfile << endl;
   }
   myfile.close();
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
      cout << states[idx][0] << ", " << states[idx][1] << ", " << 
         states[idx][2] << ", " << states[idx][1] << "\n";
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
      double x3 = (next[2] - prev[2])*((t-ti)/(tf-ti)) + prev[2];
      double x4 = (next[3] - prev[3])*((t-ti)/(tf-ti)) + prev[3];
      x.push_back(x1);
      x.push_back(x2);
      x.push_back(x3);
      x.push_back(x4);
   }
}


// applies the control from vector u to the mujoco sim
void apply_controls(mjData* d, vector<double> &u){
   d->ctrl[0] = u[0];
   d->ctrl[1] = u[1];
}

// clear all times
void cleartimers(mjData* d)
{
    for( int i=0; i<mjNTIMER; i++ )
    {
        d->timer[i].duration = 0;
        d->timer[i].number = 0;
    }
}


// show profiler
void profilerupdate(void)
{
    int i, n;

    // update constraint figure
    figconstraint.linepnt[0] = mjMIN(mjMIN(d->solver_iter, mjNSOLVER), mjMAXLINEPNT);
    for( i=1; i<5; i++ )
        figconstraint.linepnt[i] = figconstraint.linepnt[0];
    if( m->opt.solver==mjSOL_PGS )
    {
        figconstraint.linepnt[3] = 0;
        figconstraint.linepnt[4] = 0;
    }
    if( m->opt.solver==mjSOL_CG )
        figconstraint.linepnt[4] = 0;
    for( i=0; i<figconstraint.linepnt[0]; i++ )
    {
        // x
        figconstraint.linedata[0][2*i] = (float)i;
        figconstraint.linedata[1][2*i] = (float)i;
        figconstraint.linedata[2][2*i] = (float)i;
        figconstraint.linedata[3][2*i] = (float)i;
        figconstraint.linedata[4][2*i] = (float)i;

        // y
        figconstraint.linedata[0][2*i+1] = (float)d->nefc;
        figconstraint.linedata[1][2*i+1] = (float)d->solver[i].nactive;
        figconstraint.linedata[2][2*i+1] = (float)d->solver[i].nchange;
        figconstraint.linedata[3][2*i+1] = (float)d->solver[i].neval;
        figconstraint.linedata[4][2*i+1] = (float)d->solver[i].nupdate;
    }

    // update cost figure
    figcost.linepnt[0] = mjMIN(mjMIN(d->solver_iter, mjNSOLVER), mjMAXLINEPNT);
    for( i=1; i<3; i++ )
        figcost.linepnt[i] = figcost.linepnt[0];
    if( m->opt.solver==mjSOL_PGS )
    {
        figcost.linepnt[1] = 0;
        figcost.linepnt[2] = 0;
    }

    for( i=0; i<figcost.linepnt[0]; i++ )
    {
        // x
        figcost.linedata[0][2*i] = (float)i;
        figcost.linedata[1][2*i] = (float)i;
        figcost.linedata[2][2*i] = (float)i;

        // y
        figcost.linedata[0][2*i+1] = (float)mju_log10(mju_max(mjMINVAL, d->solver[i].improvement));
        figcost.linedata[1][2*i+1] = (float)mju_log10(mju_max(mjMINVAL, d->solver[i].gradient));
        figcost.linedata[2][2*i+1] = (float)mju_log10(mju_max(mjMINVAL, d->solver[i].lineslope));
    }

    // get timers: total, collision, prepare, solve, other
    int itotal = (d->timer[mjTIMER_STEP].duration > d->timer[mjTIMER_FORWARD].duration ?
                    mjTIMER_STEP : mjTIMER_FORWARD);
    float tdata[5] = { 
        (float)(d->timer[itotal].duration/mjMAX(1,d->timer[itotal].number)),
        (float)(d->timer[mjTIMER_POS_COLLISION].duration/mjMAX(1,d->timer[mjTIMER_POS_COLLISION].number)),
        (float)(d->timer[mjTIMER_POS_MAKE].duration/mjMAX(1,d->timer[mjTIMER_POS_MAKE].number)) +
            (float)(d->timer[mjTIMER_POS_PROJECT].duration/mjMAX(1,d->timer[mjTIMER_POS_PROJECT].number)),
        (float)(d->timer[mjTIMER_CONSTRAINT].duration/mjMAX(1,d->timer[mjTIMER_CONSTRAINT].number)),
        0
    };
    tdata[4] = tdata[0] - tdata[1] - tdata[2] - tdata[3];

    // update figtimer
    int pnt = mjMIN(201, figtimer.linepnt[0]+1);
    for( n=0; n<5; n++ )
    {
        // shift data
        for( i=pnt-1; i>0; i-- )
            figtimer.linedata[n][2*i+1] = figtimer.linedata[n][2*i-1];

        // assign new
        figtimer.linepnt[n] = pnt;
        figtimer.linedata[n][1] = tdata[n];
    }

    // get sizes: nv, nbody, nefc, sqrt(nnz), ncont, iter
    float sdata[6] = {
        (float)m->nv,
        (float)m->nbody,
        (float)d->nefc,
        (float)mju_sqrt((mjtNum)d->solver_nnz),
        (float)d->ncon,
        (float)d->solver_iter
    };

    // update figsize
    pnt = mjMIN(201, figsize.linepnt[0]+1);
    for( n=0; n<6; n++ )
    {
        // shift data
        for( i=pnt-1; i>0; i-- )
            figsize.linedata[n][2*i+1] = figsize.linedata[n][2*i-1];

        // assign new
        figsize.linepnt[n] = pnt;
        figsize.linedata[n][1] = sdata[n];
    }
}



// show profiler
void profilershow(mjrRect rect)
{
    mjrRect viewport = {rect.width - rect.width/5, rect.bottom, rect.width/5, rect.height/4};
    mjr_figure(viewport, &figtimer, &con);
    viewport.bottom += rect.height/4;
    mjr_figure(viewport, &figsize, &con);
    viewport.bottom += rect.height/4;
    mjr_figure(viewport, &figcost, &con);
    viewport.bottom += rect.height/4;
    mjr_figure(viewport, &figconstraint, &con);
}

//-------------------------------- simulation and rendering -----------------------------

// make option string
void makeoptionstring(const char* name, char key, char* buf)
{
    int i=0, cnt=0;

    // copy non-& characters
    while( name[i] && i<50 )
    {
        if( name[i]!='&' )
            buf[cnt++] = name[i];

        i++;
    }

    // finish
    buf[cnt] = ' ';
    buf[cnt+1] = '(';
    buf[cnt+2] = key;
    buf[cnt+3] = ')';
    buf[cnt+4] = 0;
}


// advance simulation
void simulation(void)
{
    // no model
    if( !m )
        return;

    // clear timers
    cleartimers(d);

    // paused
    if( paused )
    {
        // apply pose perturbations, run mj_forward
        if( pert.active )
        {
            mjv_applyPerturbPose(m, d, &pert, 1);      // move mocap and dynamic bodies
            mj_forward(m, d);
        }
    }

    // running
    else
    {
        // slow motion factor: 10x
        mjtNum factor = (slowmotion ? 10 : 1);

        // advance effective simulation time by 1/refreshrate
        mjtNum startsimtm = d->time;
        while( (d->time-startsimtm)*factor<1.0/refreshrate )
        {
            // clear old perturbations, apply new
            mju_zero(d->xfrc_applied, 6*m->nbody);
            if( pert.select>0 )
            {
                mjv_applyPerturbPose(m, d, &pert, 0);  // move mocap bodies only
                mjv_applyPerturbForce(m, d, &pert);
            }

            // run mj_step and count
            mj_step(m, d);

            // break on reset
            if( d->time<startsimtm )
                break;
        }
    }
}

// render
void render(GLFWwindow* window)
{
    // past data for FPS calculation
    static double lastrendertm = 0;

    // get current framebuffer rectangle
    mjrRect rect = {0, 0, 0, 0};
    glfwGetFramebufferSize(window, &rect.width, &rect.height);
    mjrRect smallrect = rect;

    // reduce rectangle when profiler is on
    if( showprofiler )
        smallrect.width = rect.width - rect.width/5;

    // no model: empty screen
    if( !m )
    {

        mjr_rectangle(rect, 0.2f, 0.3f, 0.4f, 1);
        mjr_overlay(mjFONT_NORMAL, mjGRID_TOPLEFT, rect, "Drag-and-drop model file here", 0, &con);

        // swap buffers
        glfwSwapBuffers(window); 
        return;
    }

    // advance simulation
    simulation();

    // update simulation statistics
    if( !paused )
    {
        // camera string
        char camstr[20];
        if( cam.type==mjCAMERA_FREE )
            strcpy(camstr, "Free");
        else if( cam.type==mjCAMERA_TRACKING )
            strcpy(camstr, "Tracking");
        else
            sprintf(camstr, "Fixed %d", cam.fixedcamid);

        // keyreset string
        char keyresetstr[20];
        if( keyreset<0 )
            strcpy(keyresetstr, "qpos0");
        else 
            sprintf(keyresetstr, "Key %d", keyreset);

        // solver error
        mjtNum solerr = 0;
        if( d->solver_iter )
        {
            int ind = mjMIN(d->solver_iter-1,mjNSOLVER-1);
            solerr = mju_min(d->solver[ind].improvement, d->solver[ind].gradient);
            if( solerr==0 )
                solerr = mju_max(d->solver[ind].improvement, d->solver[ind].gradient);
        }
        solerr = mju_log10(mju_max(mjMINVAL, solerr));

        // status
        sprintf(status, "%-20.1f\n%d  (%d con)\n%.3f\n%.0f\n%.2f\n%.1f  (%d it)\n%.1f %.1f\n%s\n%s\n%s\n%s",
                d->time, 
                d->nefc, 
                d->ncon,
                d->timer[mjTIMER_STEP].duration / mjMAX(1, d->timer[mjTIMER_STEP].number),
                1.0/(glfwGetTime()-lastrendertm),
                d->energy[0]+d->energy[1],
                solerr,
                d->solver_iter, 
                mju_log10(mju_max(mjMINVAL,d->solver_fwdinv[0])),
                mju_log10(mju_max(mjMINVAL,d->solver_fwdinv[1])),
                camstr, 
                mjFRAMESTRING[vopt.frame], 
                mjLABELSTRING[vopt.label],
                keyresetstr
            );
    }

    // FPS timing satistics
    lastrendertm = glfwGetTime();

    // update scene
    mjv_updateScene(m, d, &vopt, &pert, &cam, mjCAT_ALL, &scn);

    // render
    mjr_render(rect, &scn, &con);

    // show depth map
    if( showdepth )
    {
        // get the depth buffer
        mjr_readPixels(NULL, depth_buffer, rect, &con);

        // convert to RGB, subsample by 4
        for( int r=0; r<rect.height; r+=4 )
            for( int c=0; c<rect.width; c+=4 )
            {
                // get subsampled address
                int adr = (r/4)*(rect.width/4) + c/4;

                // assign rgb
                depth_rgb[3*adr] = depth_rgb[3*adr+1] = depth_rgb[3*adr+2] = 
                    (unsigned char)((1.0f-depth_buffer[r*rect.width+c])*255.0f);
            }

        // show in bottom-right corner, offset for profiler and sensor
        mjrRect bottomright = {
            smallrect.left+(3*smallrect.width)/4, 
            smallrect.bottom, 
            smallrect.width/4, 
            smallrect.height/4
        };
        mjr_drawPixels(depth_rgb, NULL, bottomright, &con);
    }

    // show overlays
    if( showhelp==1 )
        mjr_overlay(mjFONT_NORMAL, mjGRID_TOPLEFT, smallrect, "Help  ", "F1  ", &con);
    else if( showhelp==2 )
        mjr_overlay(mjFONT_NORMAL, mjGRID_TOPLEFT, smallrect, help_title, help_content, &con);

    // show info
    if( showinfo )
    {
        if( paused )
            mjr_overlay(mjFONT_NORMAL, mjGRID_BOTTOMLEFT, smallrect, "PAUSED", 0, &con);
        else
            mjr_overlay(mjFONT_NORMAL, mjGRID_BOTTOMLEFT, smallrect, 
                "Time\nSize\nCPU\nFPS\nEnergy\nSolver\nFwdInv\nCamera\nFrame\nLabel\nReset", status, &con);
    }

    // show options
    if( showoption )
    {
        int i;
        char buf[100];

        // fill titles on first pass
        if( !opt_title[0] )
        {
            for( i=0; i<mjNRNDFLAG; i++)
            {
                makeoptionstring(mjRNDSTRING[i][0], mjRNDSTRING[i][2][0], buf);
                strcat(opt_title, buf);
                strcat(opt_title, "\n");
            }
            for( i=0; i<mjNVISFLAG; i++)
            {
                makeoptionstring(mjVISSTRING[i][0], mjVISSTRING[i][2][0], buf);
                strcat(opt_title, buf);
                if( i<mjNVISFLAG-1 )
                    strcat(opt_title, "\n");
            }
        }

        // fill content
        opt_content[0] = 0;
        for( i=0; i<mjNRNDFLAG; i++)
        {
            strcat(opt_content, scn.flags[i] ? " + " : "   ");
            strcat(opt_content, "\n");
        }
        for( i=0; i<mjNVISFLAG; i++)
        {
            strcat(opt_content, vopt.flags[i] ? " + " : "   ");
            if( i<mjNVISFLAG-1 )
                strcat(opt_content, "\n");
        }

        // show
        mjr_overlay(mjFONT_NORMAL, mjGRID_TOPRIGHT, smallrect, opt_title, opt_content, &con);
    }

    // show profiler
    if( showprofiler )
    {
        if( !paused )
            profilerupdate();
        profilershow(rect);
    }

    // swap buffers
    glfwSwapBuffers(window);
}


int main(void)
{

   string MUJOCO_PATH = "/home/abajcsy/Documents/mujoco/mjpro150";
   string license = MUJOCO_PATH + "/bin/mjkey.txt";
   string times_path = "times.csv";
   string controls_path = "controls.csv";
   string states_path = "states.csv";
   string states_measured_filename = "states_measured.csv";

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

   // init GLFW
   if (!glfwInit())
      return 1;

   // get refreshrate
   refreshrate = glfwGetVideoMode(glfwGetPrimaryMonitor())->refreshRate;

   // multisampling
   glfwWindowHint(GLFW_SAMPLES, 4);

   GLFWwindow* window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
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

   vector<double> u, x;
   vector<vector<double> > states_measured;
   mjtNum start_t = d->time;
   cout << "start time: " << start_t << endl; 

   // install control callback
   //mjcb_control = mycontroller;

   // run main loop, target real-time simulation and 60 fps rendering
   while( !glfwWindowShouldClose(window))
   {
      mjtNum simstart = d->time;
      // advance interactive simulation for 1/60 sec
      //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
      //  this loop will finish on time for the next frame to be rendered at 60 fps.
      //  Otherwise add a cpu timer and exit this loop when it is time to render.
      while( d->time - simstart < 1.0/60.0 ){
         //mj_step1(m, d);
         mjtNum curr_t = d->time - start_t;

         interpolate_controls(curr_t, times, controls, u);
         apply_controls(d, u);
         cout << "control: [" << d->ctrl[0] << ", " << d->ctrl[1] << "]\n";

         //interpolate_states(curr_t, times, states, x);
         //d->qpos[0] = x[0];
         //d->qpos[1] = x[1];
         //d->qvel[0] = x[2];
         //d->qvel[1] = x[3];
         //cout << "state: [" << d->qpos[0] << ", " << d->qpos[1] << ", " << 
         //   d->qvel[0] << ", " << d->qvel[1] << "]\n";

         mj_step(m, d);

         cout << "state: [" << d->qpos[0] << ", " << d->qpos[1] << "]\n";
         vector<double> x_measured = {d->qpos[0], d->qpos[1], d->qvel[0], d->qvel[1]};
         states_measured.push_back(x_measured);
      }

      // update scene and render
      mjv_updateScene(m, d, &opt, &pert, &cam, mjCAT_ALL, &scn);
      mjr_render(viewport, &scn, &con);

      // simulate and render
      //render(window);

      // swap OpenGL buffers (blocking call due to v-sync)
      glfwSwapBuffers(window);

      // process pending GUI events, call GLFW callbacks
      glfwPollEvents();
   }

   // save out the tracked states
   write_states(states_measured_filename, states_measured);

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
