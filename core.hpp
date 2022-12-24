//
//  core.hpp
//  core
//
//  Created by Phil Huffman on 12/21/21.
//

#ifndef core_hpp
#define core_hpp

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <list>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <vector>

#define VALUE_ERROR 4
#define VECTOR_EMPTY 3 
#define DMAP_EMPTY 2 
#define GMAP_EMPTY 1 
#define OK 0


typedef unsigned long ul;
typedef unsigned long long ull;
typedef std::vector<int> vi;
typedef std::vector<double> vd;
typedef std::vector<ul> vul;
typedef std::vector<long> vl;
typedef std::map<std::string, vi> msvi;
typedef std::vector<std::string> vs;

using namespace std::chrono;

int main (int, char **);
std::string formatMicroSeconds(const ul , int , bool verbose=false );
std::string formatTime (bool doDate=false, bool doTime=true);
void init();
void randomFill(ul, vi &, std::string);
void sys(vul &, ul);

const bool FULL_Run(true);
const bool WARN_Lagards(true);
const int FORMATTED_MicroSecondLength(10);
const int GAPPER_Length(29);
const int MAX_DistroLines(150);
const int MAX_Passes(1);
const int MAX_Warnings(4);
const int MEDIAN_TrialSize(5);
const int MICROSECOND_Length(13);
const int iMax(std::numeric_limits<int>::max());
const int iMin(std::numeric_limits<int>::min());
const std::string FN_Base("/Users/prh/Keepers/code/xCode/shells/results/");
const ul MAX_SampleSize(1000000000);
const ul MIN_SampleSize(100000);
const vs DISTRO_NAMES({
                        "Bernoulli",
                        "Binomial",
                        "Gamma",
                        "Geometric",
                        "Normal",
                        "Poisson",
                        "Uniform - Sorted & Reversed",
                        "Uniform - Sorted",
                        "Uniform",
                      });
const vul SIZES({
                  100000,
                  1000000,
                  10000000,
                  100000000,
                  1000000000,
                });

struct my_numpunct : std::numpunct<char> {
  std::string do_grouping() const {return "\03";}
};

struct distroData {
  ul time;
};
typedef std::map<std::string,distroData> m_s_dd;

struct runData {
  vul gaps;
};
typedef std::map<ul,runData> m_ul_rd;

struct gs {
  int warnings;
  enum errorState {
    ok = 0,
    outOfOrder = 1,
    unknown = 1 << 31,
  } status;
  std::function<void(vul &, ul)> gapFn;
  m_ul_rd results;
};
typedef std::map<std::string,gs> m_s_gs;

struct topGapper {
  std::string gapper;
  ul time;
  topGapper(ul time, std::string gapper) {
    this -> time = time;
    this -> gapper = gapper;
  }
};
typedef topGapper tg;
typedef std::vector<tg> vtg;

struct originalSample {
  vi sample;
  vtg results;
};
typedef std::map<ul,originalSample> m_ul_os;

struct distroStruct {
  m_ul_os originals;
};
typedef std::map<std::string,distroStruct> m_s_ds;

void ciura(vul &, ul);
void frank(vul &, ul);
void gonnet(vul &, ul);
void hibbard(vul &, ul);
void knuth(vul &, ul);
void papernov(vul &, ul);
void phil(vul &, ul);
void pratt(vul &, ul);
void pratt_A(vul &, ul);
void sedgewick1985(vul &, ul);
void sedgewick82(vul &, ul);
void sedgewick86(vul &, ul);
void shell(vul &, ul);
void tokuda(vul &, ul);

#endif /* core_hpp */
