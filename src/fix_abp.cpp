/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Changed from fix_nve.cpp May 14th 2016 version. Guang Shi, July 2016
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm. 
------------------------------------------------------------------------- */

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_abp.h"    
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixABP::FixABP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg <= 6)
    error->all(FLERR,"Illegal fix abp command");

  t_target = force->numeric(FLERR,arg[3]); // set temperature
  translational_friction = force->numeric(FLERR,arg[4]); 

  translational_diff = t_target/translational_friction;

  rotational_diff = 3*translational_diff; //3D, assuming sigma=1.0

  active_velocity = force->numeric(FLERR,arg[5]); 


  seed = force->inumeric(FLERR,arg[6]); //seed for random number generator. integer


  if (t_target <= 0.0) error->all(FLERR,"Fix abp temperature must be > 0.0");
  
  if (seed <= 0) error->all(FLERR,"Illegal fix abp command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixABP::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixABP::init()
{
  dtv = update->dt;  // timestep
  dtwiener= sqrt(dtv);
  rand_trans = sqrt(2*translational_diff);
  rand_rot = sqrt(2*rotational_diff);

  // randomise orientation;
  int nlocal = atom->nlocal;
  double **orientation = atom->mu;
  double norm,theta, phi;
  
 if (fabs(orientation[0][0])<1e-6 && fabs(orientation[0][1])<1e-6 && fabs(orientation[0][2])<1e-6){
   // printf("::: Initialising a random orientation\n");
   for (int i = 0; i < nlocal; i++){
    theta=random->uniform()*M_PI;
    phi=random->uniform()*M_PI*2;
  
    orientation[i][0] = cos(phi)*sin(theta);
    orientation[i][1] = sin(phi)*sin(theta);
    orientation[i][2] = cos(theta);
    norm = sqrt(orientation[i][0]*orientation[i][0]+orientation[i][1]*orientation[i][1]+orientation[i][2]*orientation[i][2]);
    orientation[i][0]/=norm;
    orientation[i][1]/=norm;
    orientation[i][2]/=norm;
  }
 }
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixABP::initial_integrate(int vflag)
{
  double dtfm;
  double randf;

  // update v, x and angles of particles in group

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;

  double **orientation = atom->mu; 
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double e[3];//orientation
  double e_theta[3];
  double e_phi[3];
  double v0,inv_gammaTrans;
  double costheta, cosphi, sintheta,sinphi,norm,phi;
  v0 = active_velocity;
  inv_gammaTrans = 1./translational_friction;

  double r1,r2; //random numbers

  double oldx,oldy,oldz;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        e[0] = orientation[i][0];
        e[1] = orientation[i][1];
        e[2] = orientation[i][2];

        costheta = e[2];//costheta
        sintheta = sqrt(1-e[2]*e[2]);
        phi= atan2(e[1], e[0]) ;//cosphi
        cosphi= cos(phi);
        sinphi= sin(phi);

        // following Winkler et al.  Soft Matter, 2015, 11, 6680

        // unit vectors
        e_theta[0] = cosphi*costheta;
        e_theta[1] = sinphi*costheta;
        e_theta[2] = -sintheta;

        e_phi[0] = -sinphi;
        e_phi[1] = cosphi;
        e_phi[2] = 0.0;


        // copy initial positions
        oldx=x[i][0];
        oldy=x[i][1];
        oldz=x[i][2];
        // translational motion
        x[i][0] += dtv*v0*e[0]
                 + dtv*inv_gammaTrans*f[i][0]
                 + dtwiener*rand_trans*random->gaussian();
        x[i][1] += dtv*v0*e[1]
                 + dtv*inv_gammaTrans*f[i][1]
                 + dtwiener*rand_trans*random->gaussian();
        x[i][2] += dtv*v0*e[2]
                 + dtv*inv_gammaTrans*f[i][2]
                 + dtwiener*rand_trans*random->gaussian();

        // orientational motion
        r1 = random->gaussian();
        r2 = random->gaussian();
      
        e[0] += e_theta[0]*dtwiener*r1*rand_rot
              + e_phi[0]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[0];

        e[1] += e_theta[1]*dtwiener*r1*rand_rot
              + e_phi[1]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[1];
        e[2] += e_theta[2]*dtwiener*r1*rand_rot
              + e_phi[2]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[2];

        // renormalise
        norm = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
        e[0]/=norm;
        e[1]/=norm;
        e[2]/=norm;
        
        orientation[i][0]=e[0];
        orientation[i][1]=e[1];
        orientation[i][2]=e[2];

        // get the displacement
        v[i][0]=x[i][0]-oldx;
        v[i][1]=x[i][1]-oldy;
        v[i][2]=x[i][2]-oldz;
        // divide by the unit of time
        v[i][0]/=dtv;
        v[i][1]/=dtv;
        v[i][2]/=dtv;
      
      // printf("if (rmass) : to be implemented.\n");
      // exit(1);
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        e[0] = orientation[i][0];
        e[1] = orientation[i][1];
        e[2] = orientation[i][2];

        costheta = e[2];//costheta
        sintheta = sqrt(1-e[2]*e[2]);
        phi= atan2(e[1], e[0]) ;//cosphi
        cosphi= cos(phi);
        sinphi= sin(phi);

        // following Winkler et al.  Soft Matter, 2015, 11, 6680

        // unit vectors
        e_theta[0] = cosphi*costheta;
        e_theta[1] = sinphi*costheta;
        e_theta[2] = -sintheta;

        e_phi[0] = -sinphi;
        e_phi[1] = cosphi;
        e_phi[2] = 0.0;
        
        // copy initial positions
        oldx=x[i][0];
        oldy=x[i][1];
        oldz=x[i][2];

        // translational motion
        x[i][0] += dtv*v0*e[0]
                 + dtv*inv_gammaTrans*f[i][0]
                 + dtwiener*rand_trans*random->gaussian();
        x[i][1] += dtv*v0*e[1]
                 + dtv*inv_gammaTrans*f[i][1]
                 + dtwiener*rand_trans*random->gaussian();
        x[i][2] += dtv*v0*e[2]
                 + dtv*inv_gammaTrans*f[i][2]
                 + dtwiener*rand_trans*random->gaussian();

        // orientational motion
        r1 = random->gaussian();
        r2 = random->gaussian();
      
        e[0] += e_theta[0]*dtwiener*r1*rand_rot
              + e_phi[0]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[0];

        e[1] += e_theta[1]*dtwiener*r1*rand_rot
              + e_phi[1]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[1];
        e[2] += e_theta[2]*dtwiener*r1*rand_rot
              + e_phi[2]*dtwiener*r2*rand_rot
              - 2.*rotational_diff*dtv*e[2];

        // renormalise
        norm = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
        e[0]/=norm;
        e[1]/=norm;
        e[2]/=norm;
        
        orientation[i][0]=e[0];
        orientation[i][1]=e[1];
        orientation[i][2]=e[2];
             
        // get the displacement
        v[i][0]=x[i][0]-oldx;
        v[i][1]=x[i][1]-oldy;
        v[i][2]=x[i][2]-oldz;
        // divide by the unit of time
        v[i][0]/=dtv;
        v[i][1]/=dtv;
        v[i][2]/=dtv;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixABP::reset_dt()
{
  dtv = update->dt;  // timestep
  dtwiener= sqrt(dtv);
  rand_trans = sqrt(2*translational_diff);
  rand_rot = sqrt(2*rotational_diff);
}
