#include <iostream>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#define ASSERT(x) assert(x)

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "particle_tree.h"
#include "particle.h"

Particle ParticleShift(const Particle &_p, const Particle &_p_ref) {
    Particle p_new;
    p_new.mass = _p.mass;
    p_new.pos[0] = _p.pos[0] + _p_ref.pos[0];
    p_new.pos[1] = _p.pos[1] + _p_ref.pos[1];
    p_new.pos[2] = _p.pos[2] + _p_ref.pos[2];
    p_new.vel[0] = _p.vel[0] + _p_ref.vel[0];
    p_new.vel[1] = _p.vel[1] + _p_ref.vel[1];
    p_new.vel[2] = _p.vel[2] + _p_ref.vel[2];
    return p_new;
}

int main(int argc, char **argv){
    int num=0;
    int WIDTH=WRITE_WIDTH;
    int PRECISION=WRITE_PRECISION;
    int unit=0;

    int copt;
    while ((copt = getopt(argc, argv, "n:w:p:u:h")) != -1)
        switch (copt) {
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
            std::cout<<"keplertree [option] datafilename\n"
                     <<"Input data file format: hiarch_order branch_id mass1,mass2,semi,ecc,angle[3],ecc_anomaly\n"
                     <<"Example for hiarch_order & branch_id:\n"
                     <<"     -------------------------------------------------\n"
                     <<"       hiarch order          branch id                \n"
                     <<"           0                      0                   \n"
                     <<"                                 / \\                  \n"
                     <<"           1                    0   1                 \n"
                     <<"                               / \\ / \\                \n"
                     <<"           2                  0  1 2  3               \n"
                     <<"     -------------------------------------------------\n"
                     <<"     PS: if 1-0 has no children, 1-1 still have 2, 3 as children's branch_id\n"
                     <<"Options: (*) show defaulted values\n"
                     <<"   -n [int]:  number of pairs to read (defaulted: all)\n"
                     <<"   -w [int]:  print width("<<WIDTH<<")\n"
                     <<"   -p [int]:  print precision("<<PRECISION<<")\n"
                     <<"   -u [int]:  0: unscale; 1: x[PC], v[km/s], semi[AU], period[days]; 2: x[AU], v[km/s], semi[AU], period[days] (0)\n"
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

    ParticleTree<Particle> ptree;
    COMM::Binary bin;

    // unit convert 
    Float G=1.0;
    Float twopi= 8.0*std::atan(1.0);
    Float pc2au = 206264.806;
    Float kms2auyr = 0.210945021;
    if (unit>0) G = twopi*twopi; // AU^3 yr^-2 M_sun^-1

    int N=0;
    while (!fs.eof()) {
        int id,ib;
        fs>>id>>ib>>bin.m1>>bin.m2>>bin.semi>>bin.ecc>>bin.incline>>bin.rot_horizon>>bin.rot_self>>bin.ecca;
        if (fs.eof()) break;
        N++;

        if (unit==1) bin.semi *= pc2au;

        Particle p[2];
        bin.calcParticles(p[0],p[1],G);

        for (int k=0; k<2; k++) {
            for (int j=0; j<3; j++) {
                if (unit>0) p[k].vel[j] /= kms2auyr;
                if (unit==1) p[k].pos[j] /= pc2au;
            }
        }
        
        bool flag=ptree.link(id,ib,p[0],p[1],ParticleShift);
        if (!flag) {
            std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
            abort();
        }
        if (N==num) break;
    }

    N++;
    Particle plist[N];

    int count=ptree.collectAndStore(plist,N);
    if (count<0) {
        std::cerr<<"Error: particle number mismatched particle tree!\n";
        abort();
    }

    std::cout<<std::setprecision(PRECISION);
    for (int i=0; i<N; i++) {
        plist[i].printColumn(std::cout, WIDTH);
        std::cout<<std::endl;
    }
    
    return 0;
}
