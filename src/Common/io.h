#pragma once
#include <vector>

namespace COMM {

    //! IO Params container
    class IOParamsContainer{
        std::vector<double*> d_f64;
        std::vector<long long int*> d_s64;
        std::vector<int*> d_s32;
        std::vector<std::string*> d_str;
    
    public:
        void store(double* _item) {
            d_f64.push_back(_item);
        }

        void store(long long int* _item) {
            d_s64.push_back(_item);
        }

        void store(int* _item) {
            d_s32.push_back(_item);
        }

        void store(std::string* _item) {
            d_str.push_back(_item);
        }

        void writeAscii(FILE *_fout) {
            for(size_t i=0; i<d_f64.size(); i++) fprintf(_fout, "%26.15e ", *d_f64[i]);
            for(size_t i=0; i<d_s64.size(); i++) fprintf(_fout, "%lld ",    *d_s64[i]);
            for(size_t i=0; i<d_s32.size(); i++) fprintf(_fout, "%d ",      *d_s32[i]);
            for(size_t i=0; i<d_str.size(); i++) fprintf(_fout, "%s ",d_str[i]->c_str());
            fprintf(_fout,"\n");
        }
    
        void readAscii(FILE *_fin) {
            size_t rcount=0;
            for(size_t i=0; i<d_f64.size(); i++) {
                rcount=fscanf(_fin, "%lf ", d_f64[i]);
                if (rcount<1) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            for(size_t i=0; i<d_s64.size(); i++) {
                rcount=fscanf(_fin, "%lld ", d_s64[i]);
                if (rcount<1) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            for(size_t i=0; i<d_s32.size(); i++) {
                rcount=fscanf(_fin, "%d ", d_s32[i]);
                if (rcount<1) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
            }
            for(size_t i=0; i<d_str.size(); i++) {
                char dtmp[1024];
                rcount=fscanf(_fin, "%s ", dtmp);
                if (rcount<1) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                    abort();
                }
                *d_str[i] = dtmp;
            }
        }

    };

    // IO Params
    template <class Type>
    struct IOParams{
        Type value;
        const char* name;
        const char* defaulted;

        IOParams(IOParamsContainer& _ioc, const Type& _value, const char* _name, const char* _defaulted=NULL): value(_value), name(_name), defaulted(_defaulted)  {
            _ioc.store(&value);
        }

        void print(std::ostream& os) const{
            os<<name<<":   "<<value<<std::endl;
        }
    };

    template <class Type>
    std::ostream& operator <<(std::ostream& os, const IOParams<Type>& par) {
        if (par.defaulted!=NULL) os<<par.name<<": "<<par.defaulted;
        else os<<par.name<<": "<<par.value;
        return os;
    }
}
