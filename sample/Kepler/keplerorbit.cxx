#include <iostream>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#define ASSERT(x) assert(x)

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "particle.h"

int main(int argc, char **argv){
   bool iflag=false;
   int num=1;
   int WIDTH=22;
   int PRECISION=15;
   int unit=0;

   int copt;
   while ((copt = getopt(argc, argv, "in:w:p:u:h")) != -1)
       switch (copt) {
       case 'i':
           iflag=true;
           break;
       case 'n':
           num=atoi(optarg);
           break;
       case 'w':
           WIDTH=atoi(optarg);
           break;
       case 'p':
           PRECISION=atoi(optarg);
           break;
       case 'u':
           unit=atoi(optarg);
           break;
       case 'h':
           std::cout<<"keplerorbit [option] datafilename\n"
                    <<"    pariticle data (two lines): \n";
           Particle::printColumnTitle(std::cout, WIDTH);
           std::cout<<"\n    orbital data (one line): (period, r, L.(x/y/z) are not used)\n";
           COMM::BinaryTree<Particle>::printColumnTitle(std::cout, WIDTH);
           std::cout<<"\nOptions: (*) show defaulted values\n"
                    <<"   -i:        read particle data, output kepler orbit data (one line)\n"
                    <<"   -n [int]:  number of pairs(1)\n"
                    <<"   -w [int]:  print width(22)\n"
                    <<"   -p [int]:  print precision(15)\n"
                    <<"   -u [int]:  0: unscale; 1: x[PC], v[km/s], semi[AU], period[days]; 2: x[AU], v[km/s], semi[AU], period[days]\n"
                    <<std::endl;
           return 0;
       default:
           std::cerr<<"Unknown argument. check '-h' for help.\n";
           abort();
       }

   if (argc==1) {
       std::cerr<<"Please provide particle data filename\n";
       abort();
   }

   // data file name
   char* filename = argv[argc-1];

   // open data file
   std::fstream fs;
   fs.open(filename,std::fstream::in);
   if(!fs.is_open()) {
       std::cerr<<"Error: Filename "<<filename<<" not found\n";
       abort();
   }

   std::cout<<std::setprecision(PRECISION);

   Particle p[2];
   COMM::BinaryTree<Particle> bin;

   // unit convert 
   Float G=1.0;
   Float twopi= 8.0*std::atan(1.0);
   Float pc2au = 206264.806;
   Float kms2auyr = 0.210945021;
   if (unit>0) G = twopi*twopi; // AU^3 yr^-2 M_sun^-1
   //if (unit>0) G = 0.00449850214; // (km/s)^2 pc M_sun^-1

   if(iflag) {
       for(int i=0; i<num; i++) {
           p[0].readAscii(fs);
           p[1].readAscii(fs);

           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }

           bin.setMembers(&p[0],&p[1]);
           bin.calcCenterOfMass();

           // pc ->AU
           if(unit==1) {
               for (int k=0; k<2; k++) {
                   for (int j=0; j<3; j++) {
                       p[k].pos[j] *= pc2au;
                   }
               }
           }

           // km/s -> AU/yr
           if(unit>0) {
               for (int k=0; k<2; k++) {
                   for (int j=0; j<3; j++) {
                       p[k].vel[j] *= kms2auyr;
                   }
               }
           }

           bin.calcOrbit(G);
           bin.printColumn(std::cout, WIDTH);
           std::cout<<std::endl;
       }
   }
   else {
       for(int i=0; i<num; i++) {
           bin.readAscii(fs);
           if(unit==1) {
               for (int j=0; j<3; j++) {
                   bin.pos[j] *= pc2au;
               }
               bin.semi *= pc2au;
           }

           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }
           bin.calcParticles(p[0], p[1], G);
           
           for (int k=0; k<2; k++) {
               // center-of-mass correction
               for (int j=0; j<3; j++) {
                   p[k].pos[j] += bin.pos[j];
                   p[k].vel[j] += bin.vel[j];
                   if (unit>0) p[k].vel[j] /= kms2auyr;
                   if (unit==1) p[k].pos[j] /= pc2au;
               }

               p[k].printColumn(std::cout, WIDTH);
               std::cout<<std::endl;
           }
       }
   }

   return 0;
}
