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
   int unit=0;
   int width = WRITE_WIDTH;
   int precision = WRITE_PRECISION;
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
           width=atoi(optarg);
           break;
       case 'p':
           precision=atoi(optarg);
           break;
       case 'u':
           unit=atoi(optarg);
           break;
       case 'h':
           std::cout<<"keplerorbit [option] datafilename\n"
                    <<"    pariticle data (two lines): \n";
           Particle::printColumnTitle(std::cout, 12);
           std::cout<<"\n    orbital data (one line): (period, r, L.(x/y/z) are not used)\n";
           COMM::BinaryTree<Particle,COMM::Binary>::printColumnTitle(std::cout, 12);
           std::cout<<"\nOptions: (*) show defaulted values\n"
                    <<"   -i:        read particle data, output kepler orbit data (one line)\n"
                    <<"   -n [int]:  number of pairs(1)\n"
                    <<"   -w [int]:  print width(22)\n"
                    <<"   -p [int]:  print precision(15)\n"
                    <<"   -u [int]:  0: unscale \n"
                    <<"              1: m[Msun], r[AU], v[AU/yr], semi[AU], period[yr]\n"
                    <<"              2: m[Msun], r[AU], v[km/s],  semi[AU], period[days]\n"
                    <<"              3: m[Msun], r[PC], v[km/s],  semi[PC], period[days]\n"
                    <<"              4: m[Msun], r[PC], v[PC/Myr],semi[PC], period[Myr]\n";
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

   std::cout<<std::setprecision(precision);

   Particle p[2];
   COMM::BinaryTree<Particle,COMM::Binary> bin;
   bin.setMembers(&p[0],&p[1],1,2);

   // unit convert 
   Float G=1.0;
   Float twopi= 8.0*std::atan(1.0);
   Float pc2au = 206264.806;
   Float kms2auyr = 0.210945021;
   if (unit>0) G = twopi*twopi; // AU^3 yr^-2 M_sun^-1
   if (unit==4) G = 0.00449830997959438; // (pc/myr)^2 pc M_sun^-1

   if(iflag) {
       for(int i=0; i<num; i++) {
           p[0].readAscii(fs);
           p[1].readAscii(fs);

           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }

           bin.calcCenterOfMass();

           if(unit==2||unit==3) {
               // km/s -> AU/yr
               for (int k=0; k<2; k++) {
                   for (int j=0; j<3; j++) {
                       p[k].vel[j] *= kms2auyr;
                   }
               }
           }

           if(unit==3) {
               // pc ->AU
               for (int k=0; k<2; k++) {
                   for (int j=0; j<3; j++) {
                       p[k].pos[j] *= pc2au;
                   }
               }
           }

           bin.calcOrbit(G);
           // yr -> days
           if (unit>1) bin.period *= 365.25;
           bin.printColumn(std::cout, width);
           std::cout<<std::endl;
       }
   }
   else {
       for(int i=0; i<num; i++) {
           bin.readAscii(fs);
           if(unit==3) {
               // pc ->AU
               for (int j=0; j<3; j++) {
                   bin.pos[j] *= pc2au;
               }
               bin.semi *= pc2au;
           }

           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }
           bin.calcParticles(G);
           
           for (int k=0; k<2; k++) {
               // center-of-mass correction
               for (int j=0; j<3; j++) {
                   p[k].pos[j] += bin.pos[j];
                   p[k].vel[j] += bin.vel[j];
                   if (unit==2||unit==3) p[k].vel[j] /= kms2auyr;
                   if (unit==3) p[k].pos[j] /= pc2au;
               }

               p[k].printColumn(std::cout, width);
               std::cout<<std::endl;
           }
       }
   }

   return 0;
}
