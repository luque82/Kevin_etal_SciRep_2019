#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

const int P = 100000000;
const int B = 1;
const int nsteps = 750;
double PIref = 3.141592653589793238;

void boundarys(double x_p_new, double y_p_new, double z_p_new, double x_p_old, double y_p_old, double z_p_old,
               double x_left, double x_right, double y_left, double y_right, double z_upper, double z_lower, double R,
               double* x_new1, double* y_new1, double* z_new1, double* x_old1, double* y_old1, double* z_old1);

int main(int argc, char *argv[]){

//Parameters----------------------------
double x_old, x_new, y_old, y_new, z_old, z_new, wt, wx, wy, wz, Ux, Uy, Uz, U1, D0, nu, tmin, Tr, theta, alpha, tau, vx, vy, vz, phage, bac, dist, omega, Nenc, Nact ;
double v=0, Rb=1, Rp=0.1, dt=0.0066, Time=0, LX = 100, M = 1;
double x_left=-LX/2, x_right=LX/2, y_left=-LX/2, y_right=LX/2, z_lower=-LX/2, z_upper=LX/2, V = LX*LX*LX, Rcrit=Rb+Rp;

int i, p, b, n, m, k, Mf, proc, tag=314, master=0;                                                                                                                      //No. of bacteria and phage
double MatP[P][7], MatB[B][10], my_kernel[nsteps], kernel_old[nsteps], kernel_new[nsteps], en_indx[P];    //array of particle positions
MPI_Status status;

// Initialize the MPI environment
MPI_Init(&argc, &argv);

// Get the number of processes
int comm_sz;
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

// Get the rank of the process
int my_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

// Seed random numbers to rank
srand(time(NULL) + my_rank);

//Arguments
omega = strtol(argv[2], NULL, 10);   //Number of Omega's
D0 = 2.54;                           //Diffusion Constant
nu = .99;                           //Diffusion Exponent
tmin = 0.001;                        //power-law gamma parameter

Mf = M/comm_sz;                      //each worker does a chunk of iterations

printf("omega is: %f, D0 is: %f,  nu is : %f \n", omega, D0, nu);

//***************************************************************************
//                            Main Loop
//***************************************************************************

// open save locations
//FILE *pha_x = fopen("Pha_x.dat", "wb");
//FILE *pha_y = fopen("Pha_y.dat", "wb");
//FILE *pha_z = fopen("Pha_z.dat", "wb");

for(m=0;m<Mf;m++){ //iterations of experiment

//initialize old kernel
for(n=0;n<nsteps;n++){
    kernel_old[n] = 0;
}
//initialize phage activation index.
for(p=0;p<P;p++){
    en_indx[p]=1;
}

Time=0;

//***************************************************************************
//                            Initialize System
//***************************************************************************
//initialize phage
for(p=0;p<P;p++){

//random locations
double x_old=(x_right-x_left-2*Rp)*((double)rand()/(RAND_MAX)) + x_left + Rp;
double y_old=(y_right-y_left-2*Rp)*((double)rand()/(RAND_MAX)) + y_left + Rp;
double z_old=(z_upper-z_lower-2*Rp)*((double)rand()/(RAND_MAX)) + z_lower + Rp;

//store the phage tracks
MatP[p][1] = x_old; MatP[p][2] = y_old; MatP[p][3] = z_old;  //pos-new
MatP[p][4] = x_old; MatP[p][5] = y_old; MatP[p][6] = z_old;  //pos-old

//transition timer
MatP[p][7]=0;
}

//initialize bacteria
for(b=0;b<B;b++){

//random locations
double x_old=(x_right-x_left-2*Rb)*((double)rand()/(RAND_MAX)) + x_left + Rb;
double y_old=(y_right-y_left-2*Rb)*((double)rand()/(RAND_MAX)) + y_left + Rb;
double z_old=(z_upper-z_lower-2*Rb)*((double)rand()/(RAND_MAX)) + z_lower + Rb;

//store the bacteria tracks
MatB[b][1] = 0*x_old; MatB[b][2] = 0*y_old; MatB[b][3] = 0*z_old;  //pos-new
MatB[b][4] = 0*x_old; MatB[b][5] = 0*y_old; MatB[b][6] = 0*z_old;  //pos-old

//store the bacteria velocities
MatB[b][7] = 0; MatB[b][8] = 0; MatB[b][9] = 0;  //vel-new

//transition timer
MatB[b][10]=0;
}

//***************************************************************************
//                     Simulate system over nsteps
//***************************************************************************

for(n=0;n<nsteps;n++){

//Phage Dynamics
for(p=0;p<P;p++){

    if (Time >= MatP[p][7]){

    //generate random noise
        Ux = ((double)rand()/(RAND_MAX))-0.5;
        wx = sqrt(24*D0*dt)*Ux;
        Uy = ((double)rand()/(RAND_MAX))-0.5;
        wy = sqrt(24*D0*dt)*Uy;
        Uz = ((double)rand()/(RAND_MAX))-0.5;
        wz = sqrt(24*D0*dt)*Uz;

        x_old = MatP[p][4];  y_old = MatP[p][5]; z_old = MatP[p][6];

    //Langevin x-dir
        x_new = x_old + wx;
        x_old = x_new;
    //Langevin y-dir
        y_new = y_old + wy;
        y_old = y_new;
    //Langevin z-dir
        z_new = z_old + wz;
        z_old = z_new;

    //apply B.C's
    boundarys(x_new, y_new, z_new, x_old, y_old, z_old, x_left, x_right, y_left, y_right, z_upper, z_lower, Rp, &x_new, &y_new, &z_new, &x_old, &y_old, &z_old);

    //update the phage tracks
    MatP[p][1] = x_new; MatP[p][2] = y_new; MatP[p][3] = z_new;  //pos-new
    MatP[p][4] = x_old; MatP[p][5] = y_old; MatP[p][6] = z_old;  //pos-old

    //generate waiting time
    U1 = (double)rand()/(RAND_MAX);
    tau = 0*tmin*pow(1/U1,(1/nu));
    MatP[p][7] = Time + tau; //store waiting time

    }//end update

// Store phage locations
//fprintf(pha_x, "%d %f\n", n, MatB[b][1]);
//fprintf(pha_y, "%d %f\n", n, MatB[b][2]);
//fprintf(pha_z, "%d %f\n", n, MatB[b][3]);

}//end phage

//Bacteria Dynamics
for(b=0;b<B;b++){

    if (Time>=MatB[b][10]){

    //generate random run angles
        U1 = (double)rand()/(RAND_MAX);
        theta = 2*PIref*U1;
        U1 = (double)rand()/(RAND_MAX);
        alpha = 2*PIref*U1;

    //generate run time
        U1 = (double)rand()/(RAND_MAX);
        tau = -log(U1)/omega;
        MatB[b][10] = Time + tau;

    //compute velocity's
        vx=v*cos(alpha)*cos(theta);
        vy=v*cos(alpha)*sin(theta);
        vz=v*sin(alpha);
    //store the bacteria velocities
        MatB[b][7] = vx; MatB[b][8] = vy; MatB[b][9] = vz;  //vel-new
    }//end update

    x_old = MatB[b][4];  y_old = MatB[b][5]; z_old = MatB[b][6];
//Langevin x-dir
    x_new = MatB[b][7]*dt + x_old;
    x_old = x_new;
//Langevin y-dir
    y_new = MatB[b][8]*dt + y_old;
    y_old = y_new;
//Langevin z-dir
    z_new = MatB[b][9]*dt + z_old;
    z_old = z_new;

//apply B.C's
    boundarys(x_new, y_new, z_new, x_old, y_old, z_old, x_left, x_right, y_left, y_right, z_upper, z_lower, Rb, &x_new, &y_new, &z_new, &x_old, &y_old, &z_old);

//update the bacteria tracks
MatB[b][1] = x_new; MatB[b][2] = y_new; MatB[b][3] = z_new;  //pos-new
MatB[b][4] = x_old; MatB[b][5] = y_old; MatB[b][6] = z_old;  //pos-old

}//end bacteria

//*********************************Encounter Function***********************************

Nenc = 0; //set number of encounters to zero.
Nact = 0; //set number of active particles to zero.

//count number of active particles.
for(p=0;p<P;p++){
    Nact = Nact + en_indx[p];
}

//loop through phage
for(p=0;p<P;p++){

    dist = sqrt(pow((MatP[p][1]-MatB[0][1]),2)+pow((MatP[p][2]-MatB[0][2]),2)+pow((MatP[p][3]-MatB[0][3]),2));

    if(dist<Rcrit){

        Nenc = Nenc + en_indx[p];   //conditional encounter index
        en_indx[p] = 0;   // deactivate
    }

   if((fabs(MatP[p][1]-MatB[0][1])>=x_right)||(fabs(MatP[p][2]-MatB[0][2])>=y_right)||(fabs(MatP[p][3]-MatB[0][3])>=z_upper)){
        en_indx[p] = 1;   //activate
    }
} //loop on phage

kernel_new[n] = (Nenc/(M*Nact))*(V/dt);  //update kernel
kernel_old[n] = kernel_new[n];
//**************************End Encounter Function*************************//

Time = Time + dt; //update time
}// end time

}//end M

//close save locations
//fclose(bac_x);
//fclose(bac_y);
//fclose(bac_z);

if(my_rank==0){ // I am the master process and will gather results from the workers//

double kernel[nsteps];

for(n=0;n<nsteps;n++){

    kernel[n] = kernel_new[n];  //receive kernel_new from workers and update the final kernel.
    kernel_old[n] = 0;
    }//time

// open save locations
FILE *en_mat = fopen("en_mat.dat", "wb");

    for (proc=1;proc<comm_sz;proc++){

    //receive worker encounters
    MPI_Recv(&kernel_new,1,MPI_DOUBLE,proc,tag,MPI_COMM_WORLD,&status);

    for(n=0;n<nsteps;n++){

    kernel[n] = kernel_new[n] + kernel_old[n];  //receive kernel_new from workers and update the final kernel.
    kernel_old[n] = kernel[n];  //update old kernel.

    }//time

    printf("Master has received encounters from : %d\n", proc);
    }//end receive

    // Write to file
    for(n=0;n<nsteps;n++){
    fprintf(en_mat, "%d %f\n", n, kernel[n]);
    }

//close save locations
fclose(en_mat);
}
else{ //workers send encounters to the master

double tot_en=0;
    MPI_Send(&kernel_new,1,MPI_DOUBLE,master,tag,MPI_COMM_WORLD);

        for(n=0;n<nsteps;n++){
        tot_en = tot_en + kernel_new[n];
        }
    printf("Greetings from processor %d my encounters is : %f \n", my_rank, tot_en);

}

MPI_Finalize(); //let MPI finish up
return 0;
} // end main

void boundarys(double x_p_new, double y_p_new, double z_p_new, double x_p_old, double y_p_old, double z_p_old, double x_left, double x_right, double y_left, double y_right, double z_upper, double z_lower, double R, double* x_new_p, double* y_new_p, double* z_new_p, double* x_old_p, double* y_old_p, double* z_old_p){
// Apply random boundary conditions to bacteria
// R = radius of particle

    if(abs(x_p_new)>=x_right){
    x_p_new = x_p_new - (2*x_right)*x_p_new/fabs(x_p_new);
    x_p_old = x_p_old - (2*x_right)*x_p_old/fabs(x_p_old);
    }
    else if(abs(y_p_new)>=y_right){
    y_p_new = y_p_new - (2*y_right)*y_p_new/fabs(y_p_new);
    y_p_old = y_p_old - (2*y_right)*y_p_old/fabs(y_p_old);
    }
    else if(abs(z_p_new)>=z_upper){
    z_p_new = z_p_new - (2*z_upper)*z_p_new/fabs(z_p_new);
    z_p_old = z_p_old - (2*z_upper)*z_p_old/fabs(z_p_old);
    }
    *x_new_p = x_p_new;
    *y_new_p = y_p_new;
    *z_new_p = z_p_new;

    *x_old_p = x_p_old;
    *y_old_p = y_p_old;
    *z_old_p = z_p_old;

};// boundary

