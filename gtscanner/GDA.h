#ifndef _GDD_
#define _GDD_

#include <vector>
#include <map>

typedef std::vector< std::map<int, float> > GDD;

class GDA {
  private:
     std::vector<GDD> gdds;

  public:
     void newGDD(const char* file);
     GDD readOrbitFile(const char* file);

     float computeGDA(int gdd_i1, int gdd_i2);

     std::vector<GDD> getGDDs() { return gdds; }
};

#endif
