/*
The `core.cpp` file in the Shell-aio repository is a part of a shell program
written in C++. Here's a brief summary of its key functionalities:

1. **Main Function**: The main function initializes the program and returns an
OK status.

2. **Median Function**: This is a template function that calculates and returns
the median of a vector. It checks if the vector is empty and if the vector type
is arithmetic. If these conditions are met, it sorts the vector and calculates
the median.

3. **Average Function**: This is another template function that calculates and
returns the average of a vector. It also checks if the vector is empty and if
the vector type is arithmetic.

4. **ShellSort Function**: This function sorts a vector using Shell's method
with a provided gap sequence.

5. **LogError Function**: This function logs error messages.

6. **GetGreatestLength Function**: This function returns the length of the
longest string in a vector.

7. **Trim Function**: This function trims leading and trailing spaces from a
string.

8. **ZeroOrMany Function**: This function returns an "s" if the input number is
not 1, and an empty string otherwise. It's used for proper grammar in output
messages.

9. **FormatMicroSeconds Function**: This function formats microseconds into a
more readable format. It can output in verbose or compressed format.

Please note that this is a high-level summary. For a detailed understanding, you
should review the code and comments in the `core.cpp` file. l
*/

#include "core.hpp"
#include <iostream>
#include <string>

int main(int argc, char **) {
  std::locale loc(std::cout.getloc(), new my_numpunct);
  std::cout.imbue(loc);

  init();

  return OK;
}

template <typename T> T median(std::vector<T> a) {
  // returns median of a vector
  if (a.empty()) {
    logError("Median error. Empty.");
    exit(VECTOR_EMPTY);
  }

  // clangd incorrectly flags the below as an error.
  if (!std::is_arithmetic_v<T>) {
    logError("Median error. Vector type not arithmetic.");
    exit(VALUE_ERROR);
  }

  std::sort(a.begin(), a.end());
  auto mp = a.size() / 2;
  auto it(a.begin() + mp);

  if (a.size() & 1)
    return *it;

  return ((*it + *(it - 1)) / 2);
}

template <typename T> T average(const std::vector<T> a) {
  // returns mean of a vector

  if (a.empty()) {
    logError("Median error. empty vector");
    exit(VECTOR_EMPTY);
  }

  // clangd incorrectly flags the below as an error.
  if (!std::is_arithmetic_v<T>) {
    logError("Median error. Vector type not numeric.");
    exit(VALUE_ERROR);
  }

  return std::accumulate(a.begin(), a.end(), 0) / a.size();
}

template <typename T> void shellSort(std::vector<T> &v, const vul &gaps) {
  // sorts a vector using Shell's method with provided gap sequence
  for (auto gap : gaps) {
    for (auto iti(v.begin() + gap); iti != v.end(); iti++) {
      auto tmp(*iti);
      auto itj(iti);
      for (itj = iti; itj >= v.begin() + gap && *(itj - gap) > tmp; itj -= gap)
        *itj = *(itj - gap);
      *itj = tmp;
    }
  }
}

void logError(std::string msg) {
  std::string tms("");
  formatTime(tms);
  std::cerr << tms << " \t" << msg << '\n';
}

void getGreatestLength(const vs &v, int &width) {
  // returns length of longest string in vector
  width = 0;

  for (auto s : v) {
    width = width < s.size() ? s.size() : width;
  }
}

void trim(std::string &s) {
  auto fit(s.begin());
  auto rit(s.begin() + s.size() - 1);
  while (*fit == ' ') {
    fit++;
  }
  while (*rit == ' ') {
    rit--;
  }
  std::string rv(fit, rit);
  s = rv;
}

std::string zeroOrMany(ul n) { return n == 1 ? "" : "s"; }

void formatMicroSeconds(std::string &str, const ul tms, int p, bool verbose,
                        bool compresed) {
  // makes microseconds readable
  const double kd(1000000.0);
  double seconds(tms / kd);
  ul minutes(seconds / 60);
  ul secs(seconds - minutes * 60);
  seconds -= minutes * 60;
  const ul hours(minutes / 60);
  minutes -= hours * 60;

  std::stringstream sst;
  sst.precision(p);
  sst << std::setfill(' ');
  sst << std::fixed;

  if (verbose) {
    if (hours) {
      sst << hours << " hour" << zeroOrMany(hours);
    }
    if (minutes) {
      sst << (hours ? " " : "") << minutes << " minute" << zeroOrMany(minutes);
    }
    if (seconds) {
      sst << (hours || minutes ? ", and " : "") << seconds << " second"
          << zeroOrMany(secs);
    }
    str = sst.str();
    if (compresed)
      trim(str);
  }

  if (hours)
    sst << (hours > 99 ? "  " : hours > 9 ? "   " : "    ") << hours << ":";
  else
    sst << "    ";

  sst << std::setfill('0');
  if (minutes || hours)
    sst << (hours > 0     ? ""
            : minutes > 9 ? "  "
                          : "   ")
        << std::setw(minutes > 9 || hours > 0 ? 2 : 1) << minutes << ":";

  sst << (hours && minutes && secs < 10 ? "      "
          : hours && minutes            ? "     "
                                        : "")
      << std::setw((hours > 0 || minutes > 0) && secs < 10 ? 1 : 0)
      << ((hours > 0 || minutes > 0) && secs < 10 ? "0" : "") << seconds;

  str = sst.str();
  if (compresed)
    trim(str);
}

void formatTime(std::string &tms, bool doDate, bool doTime) {
  // returns date or time for display
  time_t rawtime;
  struct tm *timeinfo;
  char buffer[256];

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  std::stringstream sst;

  if (doDate) {
    strftime(buffer, 80, "%Y-%m-%d", timeinfo);
    sst << buffer;
    if (doTime)
      sst << "T";
  }
  if (doTime) {
    strftime(buffer, 80, "%H:%M:%S", timeinfo);
    sst << buffer;
  }
  tms = sst.str();
}

void gaps2string(std::string &str, const vul &gaps,
                 const std::string &delimiter) {
  // returns a numeric vector for printing
  std::stringstream sst;
  for (auto gap : gaps) {
    sst << delimiter << gap;
  }
  str = sst.str();
}

void minFillGaps(vul &gaps, ul vSize) {
  // simple gap sequence generator
  gaps.clear();
  gaps.push_back((vSize - (vSize >> 2)) | 1);
  while (gaps.back() > 1) {
    gaps.push_back((gaps.back() >> 1) | 1);
  };
}

void inspectGaps(vul &gaps, const ul vSize) {
  // check vector for readiness for use by shellSort function
  if (gaps.empty()) {
    logError("Gap error. Empty.");
    minFillGaps(gaps, vSize);
  }

  if (gaps.front() < gaps.back())
    std::reverse(gaps.begin(), gaps.end());

  if (gaps.back() != 1)
    gaps.push_back(1);

  while (gaps.front() >= vSize)
    gaps.erase(gaps.begin());

  for (auto it(gaps.begin() + 1); it != gaps.end(); it++) {
    if (*(it - 1) <= *it) {
      logError("Gap Error. Sequence.");
      minFillGaps(gaps, vSize);
    }
  }

  gaps.shrink_to_fit();
}

void writeDistros(const m_s_ds &dMap) {
  // dMap to CSV
  if (dMap.empty()) {
    logError("dMap error. Empty.");
    exit(DMAP_EMPTY);
  }

  logError("Write distros begins.");
  long maxPossibleLines(
      dMap.begin()
          ->second.originals.begin()
          ->first); // size of the first (smallest) sample vector
  long maxDistroLines(MAX_DistroLines < maxPossibleLines ? MAX_DistroLines
                                                         : maxPossibleLines);
  std::map<std::string, vi> mit;
  for (auto d0 : dMap) {
    vi tmp(d0.second.originals.begin()->second.sample.begin(),
           d0.second.originals.begin()->second.sample.begin() + maxDistroLines);
    mit[d0.first] = tmp;
  }
  std::string tms("");
  std::string fnBase(FN_Base);
  formatTime(tms);
  fnBase += tms;
  fnBase += "-distros.csv";
  std::fstream fst;
  fst.open(fnBase, std::ios::out);
  for (long indx(-1); indx < maxDistroLines; indx++) {
    for (auto vit : mit) {
      if (indx < 0) {
        fst << ',' << vit.first; // header
      } else {
        fst << ',' << vit.second[indx];
      }
    }
    fst << '\n';
  }
  fst << std::endl;
  fst.close();
  logError("Write distros ends.");
}

void writeGaps(const m_s_gs &gMap) {
  // gMap to CSV
  logError("Write gaps begins.");
  std::fstream gst;
  std::string tms(""), gs("");
  formatTime(tms);
  std::string fnBase(FN_Base + tms);
  fnBase += "-Gaps.csv";
  gst.open(fnBase, std::ios::out);
  gst << "Sequence,Size,First,Second,etc." << '\n';
  for (auto g0 : gMap) {
    for (auto g1 : g0.second.results) {
      gaps2string(gs, g1.second.gaps, ",");
      gst << g0.first << ',' << g1.first << gs << '\n';
    }
  }
  gst << std::endl;
  gst.close();
  logError("Write gaps ends.");
}

void prepG1(m_s_gs &gMap, const vul &sizes) {
  // each gapper gets a spot for each vector
  for (auto &g0 : gMap) {
    for (auto dName : DISTRO_NAMES) {
      for (auto size : sizes) {
        g0.second.results[size].gaps.clear();
      }
    }
  }
}

void makeGapSequences(m_s_gs &gMap) {
  // use attached function to produce a sequence of gap intervals
  for (auto &g0 : gMap) {
    for (auto &g1 : g0.second.results) {
      g1.second.gaps.clear();
      g0.second.gapFn(g1.second.gaps, g1.first);
      inspectGaps(g1.second.gaps, g1.first);
    }
  }
}

void make_gMap(m_s_gs &gMap, const vul &sizes) {
  // builds gMap. a spot for sequences for each size within each gapper
  // gMap["1959 Shell"].gapFn = shell;
  // gMap["1960 Frank & Lazarus"].gapFn = frank;
  // gMap["1963 Hibbard"].gapFn = hibbard;
  // gMap["1965 Papernov & Stasevich"].gapFn = papernov;
  // gMap["1971 Pratt"].gapFn = pratt;
  gMap["1973 Knuth"].gapFn = knuth;
  gMap["1982 Sedgewick"].gapFn = sedgewick82;
  gMap["1986 Sedgewick"].gapFn = sedgewick86;
  gMap["1991 Gonnet & Baeza-Yates"].gapFn = gonnet;
  gMap["1992 Tokuda"].gapFn = tokuda;
  gMap["2001 Ciura"].gapFn = ciura;
  gMap["2022 Phil"].gapFn = phil;

  prepG1(gMap, sizes);

  for (auto &g0 : gMap) {
    g0.second.status = gs::ok;
    g0.second.warnings = 0;
    for (auto &g1 : g0.second.results) {
      g1.second.gaps.clear();
    }
  }

  makeGapSequences(gMap);
  writeGaps(gMap);
  std::string tms("");
  formatTime(tms);
  std::cerr << tms << " \tgMap built. Running " << gMap.size() << " gapper"
            << zeroOrMany(gMap.size()) << ".\n";
}

void make_dMap(m_s_ds &dMap, const vul &sizes) {
  // makes dMap. a spot for samples and sort times. this sturcture can be huge
  // ensure it is passed by reference!
  std::string tms("");
  formatTime(tms);
  std::cerr << tms << " \tBuilding dMap for " << DISTRO_NAMES.size()
            << " distribution" << zeroOrMany(DISTRO_NAMES.size())
            << " each with " << sizes.size() << " sample"
            << zeroOrMany(sizes.size()) << ","
            << " large samples using complicated distributions will take time."
            << '\n';

  int maxDistroNamelength;
  getGreatestLength(DISTRO_NAMES, maxDistroNamelength);
  for (auto dName : DISTRO_NAMES) {
    for (auto size : sizes) {
      dMap[dName].originals[size].sample.clear();
      randomFill(size, dMap[dName].originals[size].sample, dName);
    }
    std::string tms("");
    formatTime(tms);
    std::cerr << tms << " \t" << std::left << std::setw(maxDistroNamelength)
              << dName << " distribution samples built.\n";
  }
  writeDistros(dMap);
  logError("dMap built.");
}

void errorFunction(const vi &wc, const vi &cc) {
  // prints a peice of an out of sequence work copy, with corresponding check
  // copy
  int n(0), w(28);

  const auto itw(wc.begin());
  const auto itc(cc.begin());
  const auto lineLimit(wc.size() > MAX_ERROR_LINES ? MAX_ERROR_LINES
                                                   : wc.size());

  std::cerr << std::right << std::setw(4) << "n" << std::right << std::setw(w)
            << "expected" << std::right << std::setw(w) << "result" << '\n';
  while (n < lineLimit) {
    std::cerr << std::right << std::setw(3) << n << std::right << std::setw(w)
              << *(itc + n) << std::right << std::setw(w) << *(itw + n) << '\n';
    n++;
  }
}

void writeTimes(m_s_ds dMap) {
  // sort times to CSV
  if (FULL_Run) {
    std::fstream fst;
    std::string tms("");
    formatTime(tms);
    std::string fileName(FN_Base + tms);
    fileName += "-Times.csv";
    fst.open(fileName, std::ios::out);
    fst << "Distro,Size,Sequence,Time\n";
    for (auto d0 : dMap) {
      for (auto d1 : d0.second.originals) {
        sort(d1.second.results.begin(), d1.second.results.end(),
             [](tg &lhs, tg &rhs) { return lhs.time < rhs.time; });
        for (auto result : d1.second.results) {
          fst << d0.first << ',' << d1.first << ',' << result.gapper << ','
              << result.time << '\n';
        }
      }
    }
    fst << std::endl;
    fst.close();
  }
}

void summerize(m_s_ds &dMap, m_s_gs gMap) {
  // prints sort times
  for (auto &d0 : dMap) {
    for (auto &d1 : d0.second.originals) {
      auto firstTime(true);
      for (auto result : d1.second.results) {
        if (firstTime) {
          firstTime = false;
          std::cout << " \nn: " << std::right << std::setw(11) << d1.first
                    << " \tDistribution: " << d0.first << '\n';
        }
        std::string fms(""), gs("");
        formatMicroSeconds(fms, result.time, 3, false, true);
        gaps2string(gs, gMap[result.gapper].results[d1.first].gaps, " ");
        std::cout << std::right << std::setw(FORMATTED_MicroSecondLength) << fms
                  << std::right << std::setw(GAPPER_Length) << result.gapper
                  << " Gaps:" << gs << '\n';
      }
    }
    std::cout << '\n';
  }
}

void listWinners(m_s_ds &dMap) {
  // prints quickest sequence for each distro and sample size
  if (FULL_Run) {
    for (auto &d0 : dMap) {
      for (auto &d1 : d0.second.originals) {
        sort(d1.second.results.begin(), d1.second.results.end(),
             [](tg &lhs, tg &rhs) { return lhs.time < rhs.time; });
        std::cout << "Quickest gapper for distro " << d0.first
                  << " with a size of " << d1.first << " is "
                  << d1.second.results.front().gapper << " which used "
                  << d1.second.results.front().time << "µs\n";
      }
      std::cout << '\n';
    }
  }
}

void eoj(m_s_gs &gMap, m_s_ds &dMap) {
  // end of job logic goes here
  writeTimes(dMap);
  summerize(dMap, gMap);
  listWinners(dMap);
}

void doSort(std::pair<const std::string, distroStruct> &d0,
            std::pair<const unsigned long, originalSample> &d1,
            std::pair<const std::string, gs> &g0, vul &gTimes, const ul size,
            vi checkCopy) {
  // this is where shellSort is called
  for (auto &g1 : g0.second.results) {
    if (g1.first == d1.first) {
      vl times;
      times.clear();
      while (times.size() < MEDIAN_TrialSize) {
        auto workCopy(d1.second.sample);
        auto start = high_resolution_clock::now();
        shellSort(workCopy, g0.second.results[d1.first].gaps);
        auto stop = high_resolution_clock::now();
        auto durT = duration_cast<microseconds>(stop - start).count();
        if (workCopy == checkCopy) {
          times.push_back(durT);
        } else {
          std::string tms("");
          formatTime(tms);
          std::cerr << tms << d0.first << " \t" << g0.first << " \t" << g1.first
                    << " \tSort error\n";
          g0.second.status = gs::outOfOrder;
          g0.second.warnings = iMax;
          errorFunction(workCopy, checkCopy);
          break;
        }
      }
      auto dur(times.size() == 1 ? times.front() : median(times));
      std::string fms("");
      formatMicroSeconds(fms, dur, 3, true, true);
      std::string tms("");
      formatTime(tms);
      std::cout << tms << std::right << std::setw(MICROSECOND_Length) << dur
                << "µs "
                // << std::right
                // << std::setw(FORMATTED_MicroSecondLength)
                << fms << " \t" << g0.first << '\n';
      d1.second.results.push_back(tg(dur, g0.first));
      gTimes.push_back(dur);
    }
  }
}

void checkForLagards(vtg results, m_s_gs &gMap, const vul &gTimes,
                     ul warnLimit) {
  // gigs slower gap sequences and "blesses" ones that have been giged if they
  // were faster than average
  auto kudo(average(gTimes));
  auto limt(kudo + (kudo >> 2));
  std::string tms("");
  for (auto result : results) {
    if (gMap[result.gapper].warnings < warnLimit) {
      if (result.time > limt) {
        gMap[result.gapper].warnings++;
        formatTime(tms);
        std::cerr << tms << " \tWarned " << result.gapper << " ("
                  << gMap[result.gapper].warnings << "/" << warnLimit << ")\n";
      } else if (result.time < kudo && gMap[result.gapper].warnings > 0) {
        gMap[result.gapper].warnings--;
        formatTime(tms);
        std::cerr << tms << " \tBlessed " << result.gapper << " ("
                  << gMap[result.gapper].warnings << "/" << warnLimit << ")\n";
      }
    }
  }
}

void work(m_s_gs &gMap, m_s_ds &dMap) {
  // traverses each gapper, sample, distro
  auto warnLimit(SIZES.size() - 1 < MAX_Warnings ? MAX_Warnings
                                                 : SIZES.size() - 1);
  for (auto &d0 : dMap) {                  // each distro
    for (auto &d1 : d0.second.originals) { // each sample
      std::string tms("");
      formatTime(tms);
      std::cout << "\n\n"
                << tms << " \tn: " << d1.first << " distro: " << d0.first
                << '\n';
      auto checkCopy(d1.second.sample);
      std::sort(checkCopy.begin(), checkCopy.end());
      vul gTimes;
      for (auto &g0 : gMap) {                 // each sequence
        if (g0.second.warnings < warnLimit) { // skip slow sequences
          doSort(d0, d1, g0, gTimes, d1.first, checkCopy);
        } else {
          std::string tms("");
          formatTime(tms);
          std::cerr << tms << " \tSkipping " << g0.first << " (too slow)\n";
        }
      }
      if (WARN_Lagards) {
        checkForLagards(d1.second.results, gMap, gTimes, warnLimit);
      }
      sort(d1.second.results.begin(), d1.second.results.end(),
           [](tg &lhs, tg &rhs) { return lhs.time < rhs.time; });
      formatTime(tms);
      std::cout << tms << " \tBest sequence for size of " << d1.first
                << " with distro " << d0.first << " is "
                << d1.second.results.front().gapper << '\n';
    }
  }
}

void clearResults(m_s_ds &dMap, m_s_gs &gMap) {
  for (auto &d0 : dMap) {
    for (auto &d1 : d0.second.originals) {
      d1.second.results.clear();
    }
  }
  for (auto &g0 : gMap) {
    g0.second.warnings = 0;
  }
}

void init() {
  // initialization function sets up dMap & gMap for sorting
  auto t0(system_clock::now());
  std::string tms("");
  formatTime(tms);
  std::cerr << tms << " \tInitialization begins."
            << " \tSample size for median measurement: " << MEDIAN_TrialSize
            << "\n";
  auto sizes(SIZES);
  for (auto &size : sizes) {
    size = size > MAX_SAMPLE_SIZE   ? MAX_SAMPLE_SIZE
           : size < MIN_SAMPLE_SIZE ? MIN_SAMPLE_SIZE
                                    : size;
  }

  // sizes must be unique because it is used as an index in both dMap & gMap.
  std::sort(sizes.begin(), sizes.end());
  sizes.erase(std::unique(sizes.begin(), sizes.end()),
              sizes.end()); // ensure sizes are unique

  m_s_ds dMap;
  m_s_gs gMap;

  make_gMap(gMap, sizes);
  make_dMap(dMap, sizes);
  auto t1(system_clock::now());
  auto durT = duration_cast<microseconds>(t1 - t0).count();
  std::string fms("");
  formatMicroSeconds(fms, durT, 1, true, true);
  formatTime(tms);
  std::cerr << tms << " \tInitialization complete. Required " << fms << ".\n";

  if (dMap.empty()) {
    logError("dMap uninitalized.");
    exit(DMAP_EMPTY);
  }

  if (gMap.empty()) {
    logError("gMap uninitalized.");
    exit(GMAP_EMPTY);
  }

  for (int passes(0); passes < MAX_Passes; passes++) {
    if (FULL_Run)
      work(gMap, dMap);

    auto t2(system_clock::now());
    auto durE = duration_cast<microseconds>(t2 - t1).count();
    t1 = system_clock::now();
    std::string fms(""), tms("");
    formatMicroSeconds(fms, durE, 1, true, true);
    formatTime(tms);
    std::cerr << tms << " \tPass " << (passes + 1) << "/" << MAX_Passes
              << " complete. Required " << fms << ".\n";
    clearResults(dMap, gMap);
  }
  eoj(gMap, dMap);
}

// ===== gap sequence generator functions that go into gMap ===========
void shell(vul &gaps, ul vSize) {
  ul gap(vSize);
  while (gap > 1) {
    gap >>= 1;
    gaps.push_back(gap);
  }
}

void frank(vul &gaps, ul vSize) {
  ul gap(vSize >> 1);
  while (gap) {
    gaps.push_back(gap | 1);
    gap >>= 1;
  }
}

void hibbard(vul &gaps, ul vSize) {
  ul gap(1);
  while (gap < vSize) {
    gaps.push_back(gap);
    gap <<= 1;
    gap |= 1;
  }
}

void papernov(vul &gaps, ul vSize) {
  ul n(1);
  gaps.push_back(n);
  while (gaps.back() < vSize)
    gaps.push_back((2 << n++) + 1);
  gaps.pop_back();
}

bool is3smooth(ul n) {
  while (n % 3 == 0)
    n /= 3;
  while (n % 2 == 0)
    n /= 2;
  return n == 1;
}

void pratt(vul &gaps, ul vSize) {
  gaps.clear();
  for (ul n(1); n < vSize; n++)
    if (is3smooth(n))
      gaps.push_back(n);
}

void knuth(vul &gaps, ul vSize) {
  ul k(1), lim(vSize / 3);
  do {
    k *= 3;
    gaps.push_back((k - 1) >> 1);
  } while (gaps.back() < lim);
  gaps.pop_back();
}

bool mySeq(ul a, ul b) { return a > b; }

ul pw2(ul e) {
  ul rv(1);
  while (e--) {
    rv *= 2;
  }
  return rv;
}

void sedgewick82(vul &gaps, ul vSize) {
  gaps.push_back(1);
  for (ul k(1); gaps.back() < vSize; k++) {
    gaps.push_back(pw2(k + 1) + 3 * pw2(k - 1) + 1);
  }
}

void sedgewick86(vul &gaps, ul vSize) {
  ul k(0);
  gaps.push_back(1);
  do {
    if (++k & 1) {
      gaps.push_back((2 << (k + 3)) - (6 * (1 << (k >> 1)) + 1));
    } else {
      gaps.push_back(9 * (2 << k) - (1 << k) + 1);
    }
  } while (gaps.back() < vSize);
  gaps.pop_back();
}

void gonnet(vul &gaps, ul vSize) {
  ul k(vSize);
  while (k > 1) {
    k = (5 * k - 1) / 11;
    k = k > 1 ? k : 1;
    gaps.push_back(k);
  }
}

void tokuda(vul &gaps, ul vSize) {
  gaps.push_back(1);
  for (ul i(1); gaps.back() < vSize; i++) {
    double a(pow(2.25, static_cast<double>(i)));
    ul j((9.0 * a - 4.0) / 5.0);
    gaps.push_back(j | 1);
  }
}

void ciura(vul &gaps, ul vSize) {
  vul t{1, 4, 10, 23, 57, 132, 301, 701, 1750};
  gaps = t;
  while (gaps.back() < vSize / 3) {
    gaps.push_back(gaps.back() << 3);
  }
}

void gap22(vul &gaps, ul vSize, int bits, vi sri, vi srj) {
  // generates gap sequence based on parameters from function phil
  bits = bits < 1 ? 1 : bits > 6 ? 6 : bits;
  auto ladd(1);
  while (--bits > 0) {
    ladd <<= 1;
    ladd |= 1;
  }
  long tmp(vSize);
  for (auto i : sri) {
    tmp -= vSize >> i;
  }
  gaps.push_back(tmp < vSize ? tmp > 0 ? tmp | 1 : vSize >> 2 : vSize >> 1);
  while (gaps.back() > ladd) {
    long tmp(gaps.back());
    for (auto j : srj) {
      long t0(gaps.back() >> j);
      if (t0 <= ladd)
        break;
      tmp -= gaps.back() >> j;
    }
    if (tmp >= gaps.back())
      break;
    gaps.push_back(tmp | ladd);
  }
  if (gaps.back() != 1)
    gaps.push_back(1);
}

void phil(vul &gaps, ul vSize) {
  int bits(4); // size of mask
  vi sri({3, 4, 5});
  vi srj({1, 3, 10});
  gap22(gaps, vSize, bits, sri, srj);
}

// ===== sort sample fillers =================================================

void getRandyBe(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution dist;
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyBi(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::binomial_distribution<int> dist;
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyGe(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::geometric_distribution<int> dist;
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyG(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<double> dist;
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyN(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dist(0, 100);
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyP(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::poisson_distribution<int> dist;
  while (n--) {
    v.push_back(dist(rd));
  }
}

void getRandyU(vi &v, ul n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(iMin, iMax);
  while (n--) {
    v.push_back(dist(gen));
  }
}

void randomFill(ul size, vi &v, std::string distroName) {
  // fills sample with specified size & random distribution
  v.clear();
  if (distroName == "Bernoulli") {
    getRandyBe(v, size);
  } else if (distroName == "Binomial") {
    getRandyBi(v, size);
  } else if (distroName == "Geometric") {
    getRandyGe(v, size);
  } else if (distroName == "Gamma") {
    getRandyG(v, size);
  } else if (distroName == "Normal") {
    getRandyN(v, size);
  } else if (distroName == "Poisson") {
    getRandyP(v, size);
  } else if (distroName == "Uniform") {
    getRandyU(v, size);
  } else if (distroName == "Uniform - Sorted") {
    getRandyU(v, size);
    std::sort(v.begin(), v.end());
  } else if (distroName == "Uniform - Sorted & Reversed") {
    getRandyU(v, size);
    std::sort(v.begin(), v.end());
    std::reverse(v.begin(), v.end());
  } else {
    std::cerr << "Unknown distribution requested (" << distroName
              << "). Using uniform." << std::endl;
    getRandyU(v, size);
  }
  v.shrink_to_fit();
}
