#pragma once
#include <iostream>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

namespace COMM {

#define COMM_PRINT_WIDTH 15
#define COMM_PRINT_PRECISION 7

// Use C99 hex float text for exact double round-trip in ASCII parameter files.

//! Print format parameters for aligned help-text output
struct IOParamsPrintHelp{
    int offset_short_key;
    int offset_long_key;
    int width_key;

    IOParamsPrintHelp(const int _offset_short_key, const int _offset_long_key, const int _width_key): 
        offset_short_key(_offset_short_key), offset_long_key(_offset_long_key), width_key(_width_key) {}

    static void printTypeShortNameDescription(std::ostream& os) {
        os<<"F: 64bit floating (decimal or C99 hex-float text); "
          <<"D: 32bit integer; I: 64bit integer (long); "
          <<"L: 64bit integer (long long); S: string\n";
    }

    static char getValueTypeShortName(const int& value) {
        return 'D';
    }

    static char getValueTypeShortName(const long int& value) {
        return 'I';
    }

    static char getValueTypeShortName(const long long int& value) {
        return 'L';
    }

    static char getValueTypeShortName(const double& value) {
        return 'F';
    }

    static char getValueTypeShortName(const std::string& value) {
        return 'S';
    }
};


//! IO Params — a named parameter with value, CLI key, description, and help-printing support.
/*!
    Supports two construction styles:
    - PeTar style: (container, value, key, description, defaulted, print_help_flag)
    - SDAR 3-arg compat: (container, value, name)            → description = key = name
*/
template <class Type>
struct IOParams{
    Type value;
    const char* key;
    const char* description;
    const char* defaulted;
    bool print_help_flag;

    //! Full constructor (6 args) — PeTar style with all fields.
    template <class TContainer>
    IOParams(TContainer& _ioc, const Type& _value, const char* _key, const char* _description,
             const char* _defaulted, bool _print_help_flag)
        : value(_value), key(_key), description(_description),
          defaulted(_defaulted), print_help_flag(_print_help_flag) {
        _ioc.store(_key, this);
    }

    //! PeTar 5-arg: (container, value, key, description, defaulted)
    template <class TContainer>
    IOParams(TContainer& _ioc, const Type& _value, const char* _key, const char* _description,
             const char* _defaulted)
        : IOParams(_ioc, _value, _key, _description, _defaulted, true) {}

    //! PeTar 4-arg: (container, value, key, description)
    template <class TContainer>
    IOParams(TContainer& _ioc, const Type& _value, const char* _key, const char* _description)
        : IOParams(_ioc, _value, _key, _description, (const char*)nullptr, true) {}

    //! SDAR 3-arg compat: (container, value, name) → description = key = name
    template <class TContainer>
    IOParams(TContainer& _ioc, const Type& _value, const char* _name)
        : IOParams(_ioc, _value, _name, _name, (const char*)nullptr, true) {}

    void print(std::ostream& os) const{
        std::stringstream ss(description);
        std::string token;
        std::getline(ss, token, ';');
        os<<token<<": "<<value<<std::endl;
    }
    
    void printHelp(std::ostream& os, const IOParamsPrintHelp& _align,
                   const bool print_short_flag, const bool always_print=false) const {
        if (print_help_flag || always_print) {
            bool print_flag = false;
            int multiline_width = 0;

            if (strlen(key)==1 && print_short_flag) {
                os<<std::setw(_align.offset_short_key)<<"-"<<key<<"  ";
                multiline_width = _align.offset_short_key+5;
                print_flag = true;
            }
            else if (strlen(key)>1 && !print_short_flag) {
                os<<std::setw(_align.offset_long_key)<<"--"
                  <<std::left<<std::setw(_align.width_key)<<key<<std::right;   
                multiline_width = _align.offset_long_key+_align.width_key+5;
                print_flag = true;
            }

            if (print_flag) {
                os<<"["<<IOParamsPrintHelp::getValueTypeShortName(value)<<"] ";

                std::vector<std::string> tokens;
                std::stringstream ss(description);
                std::string token;
                while (std::getline(ss, token, ';')) {
                    tokens.push_back(token);
                }

                if (tokens.size()>1) {
                    if (defaulted!=NULL) os<<tokens[0]<<": "<<defaulted<<std::endl;
                    else os<<tokens[0]<<": "<<value<<std::endl;
                    for (size_t i=1; i<tokens.size(); i++) 
                        os<<std::setw(multiline_width)<<" "<<tokens[i]<<std::endl;
                }
                else os<<*this<<std::endl;
            }
        }
    }
};

template <class Type>
std::ostream& operator <<(std::ostream& os, const IOParams<Type>& par) {
    if (par.defaulted!=NULL) os<<par.description<<": "<<par.defaulted;    
    else os<<par.description<<": "<<par.value;
    return os;
}

//! IO Params container
/*!
    Stores IOParams by name-keyed maps for type-safe read/write.
    Uses a type-prefixed ASCII format for robust round-trip serialization.
*/
class IOParamsContainer{
    struct char_cmp {
        bool operator () (const char *a,const char *b) const {
            return strcmp(a,b)<0;
        }
    };
    std::map<const char*, IOParams<double>*, char_cmp> d_f64;
    std::map<const char*, IOParams<long int>*, char_cmp> d_l64;
    std::map<const char*, IOParams<long long int>*, char_cmp> d_ll64;
    std::map<const char*, IOParams<int>*, char_cmp> d_s32;
    std::map<const char*, IOParams<std::string>*, char_cmp> d_str;
    std::map<const char*, int, char_cmp> name_types;
    
public:
    void store(const char* _name, IOParams<double>* _item) {
        d_f64[_name] = _item;
        name_types[_name] = 0;
    }

    void store(const char* _name, IOParams<long int>* _item) {
        d_l64[_name] = _item;
        name_types[_name] = 1;
    }
    
    void store(const char* _name, IOParams<long long int>* _item) {
        d_ll64[_name] = _item;
        name_types[_name] = 2;
    }

    void store(const char* _name, IOParams<int>* _item) {
        d_s32[_name] = _item;
        name_types[_name] = 4;
    }

    void store(const char* _name, IOParams<std::string>* _item) {
        d_str[_name] = _item;
        name_types[_name] = 3;
    }

    void writeAscii(FILE *_fout) {
        for(auto iter = d_f64.begin(); iter!=d_f64.end(); iter++) {
            fprintf(_fout, "%c %s %a\n", IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        }
        for(auto iter = d_l64.begin(); iter!=d_l64.end(); iter++)
            fprintf(_fout, "%c %s %ld\n",  IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        for(auto iter = d_ll64.begin();iter!=d_ll64.end();iter++)
            fprintf(_fout, "%c %s %lld\n", IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        for(auto iter = d_s32.begin(); iter!=d_s32.end(); iter++)
            fprintf(_fout, "%c %s %d\n",   IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value);
        for(auto iter = d_str.begin(); iter!=d_str.end(); iter++)
            fprintf(_fout, "%c %s %s\n",   IOParamsPrintHelp::getValueTypeShortName(iter->second->value), iter->first, iter->second->value.c_str());
    }
    
    void readAscii(FILE *_fin) {        
        char key_name[1024];
        while (!feof(_fin)) {
            size_t rcount=0;
            char type_id;
            rcount=fscanf(_fin, "%c ", &type_id);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }

            switch (type_id) {
            case 'F':
            {
                char dtmp_str[1024];
                rcount=fscanf(_fin, "%s %1023s\n", key_name, dtmp_str);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                char* end_ptr = nullptr;
                double dtmp = std::strtod(dtmp_str, &end_ptr);
                if (end_ptr == dtmp_str || *end_ptr != '\0') {
                    std::cerr<<"Error: floating parameter value '"<<dtmp_str<<"' for key "<<key_name<<" cannot be parsed.\n";
                    abort();
                }
                auto search = d_f64.find(key_name);
                if (search == d_f64.end())
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'I':
            {
                long int dtmp;
                rcount=fscanf(_fin, "%s %ld\n", key_name, &dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_l64.find(key_name);
                if (search == d_l64.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'L':
            {
                long long int dtmp;
                rcount=fscanf(_fin, "%s %lld\n", key_name, &dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_ll64.find(key_name);
                if (search == d_ll64.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'D':
            {
                int dtmp;
                rcount=fscanf(_fin, "%s %d\n", key_name, &dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_s32.find(key_name);
                if (search == d_s32.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            case 'S':
            {
                char dtmp[1024];
                rcount=fscanf(_fin, "%s %s\n", key_name, dtmp);
                if (rcount<2) {
                    std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
                    abort();
                }
                auto search = d_str.find(key_name);
                if (search == d_str.end()) 
                    std::cerr<<"Warning: parameter name key "<<key_name<<" is not found!\n";
                else 
                    search->second->value = dtmp;
                break;
            }
            default:
                std::cerr<<"Warning: parameter type not found, given "<<type_id<<", should be one of F, D, I, L, S\n";
                break;
            }
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL        
    void mpi_broadcast() {
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++)
            MPI_Bcast(&(iter->second->value), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for(auto iter=d_l64.begin(); iter!=d_l64.end(); iter++)
            MPI_Bcast(&(iter->second->value), 1, MPI_LONG, 0, MPI_COMM_WORLD);
        for(auto iter=d_ll64.begin(); iter!=d_ll64.end(); iter++)
            MPI_Bcast(&(iter->second->value), 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        for(auto iter=d_s32.begin(); iter!=d_s32.end(); iter++)
            MPI_Bcast(&(iter->second->value), 1, MPI_INT, 0, MPI_COMM_WORLD);
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) {
            size_t str_size = iter->second->value.size();
            unsigned long long str_size_ull = str_size;
            MPI_Bcast(&str_size_ull, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            str_size = static_cast<size_t>(str_size_ull);
            iter->second->value.resize(str_size);
            if (str_size > 0) {
                MPI_Bcast(&(iter->second->value[0]), str_size, MPI_CHAR, 0, MPI_COMM_WORLD);
            }
        }
    }
#endif

    void print(std::ostream& os) const{
        for(auto iter=d_f64.begin(); iter!=d_f64.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_l64.begin(); iter!=d_l64.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_ll64.begin(); iter!=d_ll64.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_s32.begin(); iter!=d_s32.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
        for(auto iter=d_str.begin(); iter!=d_str.end(); iter++) os<<iter->first<<": "<<iter->second->value<<std::endl;
    }

    void printHelp(std::ostream& os, const bool print_format_info_flag=true,
                   const bool print_all_flag=false, const int _offset_short_key=3,
                   const int _offset_long_key=4, const int _width_key=23) const{
        IOParamsPrintHelp print_help(_offset_short_key, _offset_long_key, _width_key);
        if (print_format_info_flag) {
            os<<"** Default values are shown after ':'\n"
              <<"** The char in [] indicates argument type: ";
            print_help.printTypeShortNameDescription(os);
        }
        for (int i = 0; i < 2; i++) {
            for (auto iter=name_types.begin(); iter!=name_types.end(); iter++) {
                auto name = iter->first;
                auto type = iter->second;
                bool print_short_flag = (i==0);
                if (type == 0) {
                    auto search = d_f64.find(name);
                    if (search != d_f64.end()) search->second->printHelp(os, print_help, print_short_flag, print_all_flag);
                }
                else if (type == 1) {
                    auto search = d_l64.find(name);
                    if (search != d_l64.end()) search->second->printHelp(os, print_help, print_short_flag, print_all_flag);
                }
                else if (type == 2) {
                    auto search = d_ll64.find(name);
                    if (search != d_ll64.end()) search->second->printHelp(os, print_help, print_short_flag, print_all_flag);
                }
                else if (type == 4) {
                    auto search = d_s32.find(name);
                    if (search != d_s32.end()) search->second->printHelp(os, print_help, print_short_flag, print_all_flag);
                }
                else if (type == 3) {
                    auto search = d_str.find(name);
                    if (search != d_str.end()) search->second->printHelp(os, print_help, print_short_flag, print_all_flag);
                }
                else {
                    std::cerr<<"Warning: parameter name "<<name<<" has unknown type "<<type<<", should be one of 0, 1, 2, 3, 4\n";
                }
            }
        }
    }

    //! check whether the key is defined
    bool isDefined(const char* _key) const {
        auto search = name_types.find(_key);
        if (search != name_types.end()) return true;
        return false;
    }    
};

//! check if the options are defined
/*! If option is not defined, print error message and abort

    @param[in] io_par_list list of IOParamsContainer
    @param[in] argc number of arguments
    @param[in] argv argument list
*/
static void FindUndefinedOptions(std::vector<IOParamsContainer*> io_par_list, const int argc, char* argv[], std::vector<std::string>* known_options=NULL) {
    for (int i=1; i<argc; i++) {
        if (argv[i][0]=='-') {
            std::string arg(argv[i]);
            if (arg[0]=='-' && arg[1]=='-') {
                arg = arg.substr(2);
            }
            else if (arg[0]=='-') {
                arg = arg.substr(1);
                // exclude negative number arg with '-'
                if (arg[0]>='0' && arg[0]<='9') continue;
            }
            bool found = false;
            if (known_options != NULL) {
                for (auto iter = known_options->begin(); iter != known_options->end(); ++iter) {
                    if (arg == *iter) {
                        found = true;
                        break;
                    }
                }
            }
            for (auto iter = io_par_list.begin(); iter != io_par_list.end(); ++iter) { 
                if ((*iter)->isDefined(arg.c_str())) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr<<"Error: option "<<arg<<" is not defined!\n";
                abort();
            }
        }
    }
}

}
