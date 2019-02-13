#include <iostream>

#include "particle.h"
#include "particle_group.h"

using namespace AR;

struct PtclTest{
    int id;
    Float pad;

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fout) const {
        fwrite(this, sizeof(*this),1,_fout);
    }


    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this),1,_fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! write class data to file with ASCII format
    /*! @param[in] _fout: std:osteram file for output
     */
    void writeAscii(std::ostream& _fout) const {
        _fout<<id<<" ";
    }

    //! read class data to file with ASCII format
    /*! @param[in] _fin: std::istream file for input
     */
    void readAscii(std::istream&  _fin) {
        _fin>>id;
    }
    
    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"id";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<id;
    }
    
};

void printPtclAdr(ParticleGroup<PtclTest,PtclTest>& _group) {
    const int num = _group.getSize();
    std::cout<<"Particle original address\n";
    for (int i=0; i<num; i++) {
        std::cout<<_group.getMemberOriginAddress(i)<<" ";
    }
    std::cout<<std::endl;
}

int main(int argc, char **argv){
    
    const int num = 10;
    PtclTest ptcl[num];

    std::cout<<"Sample particles:\n";
    for (int i=0; i<num; i++) {
        ptcl[i].id = i;
        std::cout<<ptcl[i].id<<" ";
    }
    std::cout<<std::endl;


    ParticleGroup<PtclTest,PtclTest> group;

    group.setMode(ListMode::copy);
    group.reserveMem(num);
    
    std::cout<<"Test addParticle\n";
    for (int i=0; i<num; i++) {
        group.addMemberAndAddress(ptcl[i]);
    }
    group.printColumn(std::cout, 4);
    std::cout<<std::endl;

    std::cout<<"Test removeParticle, list:\n";
    int rmlist[4] = {0,7,3,4};
    for (int i=0; i<4; i++) std::cout<<rmlist[i]<<" ";
    std::cout<<std::endl;

    group.removeMemberList(rmlist,4);

    std::cout<<"After remove\n";
    group.printColumn(std::cout, 4);
    std::cout<<std::endl;
    printPtclAdr(group);

    for (int i=0; i<4; i++) group.addMemberAndAddress(ptcl[rmlist[i]]);
    std::cout<<"Add back\n";
    group.printColumn(std::cout, 4);
    std::cout<<std::endl;
    printPtclAdr(group);
    
    int rmlistall[num];
    for (int i=0; i<num; i++) rmlistall[i] = num-i-1;
    std::cout<<"Inverse order remove list:\n";
    for (int i=0; i<num; i++) std::cout<<rmlistall[i]<<" ";
    std::cout<<std::endl;
    
    group.removeMemberList(rmlistall,num);
    std::cout<<"After remove\n";
    group.printColumn(std::cout, 4);
    std::cout<<std::endl;

    group.clear();
    group.setMode(ListMode::link);
    group.linkMemberArray(ptcl, num);
    std::cout<<"Clear and link to particle\n";
    group.printColumn(std::cout, 4);
    std::cout<<std::endl;

    return 0;
}
