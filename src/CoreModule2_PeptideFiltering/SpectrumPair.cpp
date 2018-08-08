// Header Include
#include "Logger.h"
#include "SpectrumPair.h"

using namespace std;

namespace specnets
{
  SpectrumPair::SpectrumPair()
    : spec1(0), spec2(0), score1(0.0), score2(0.0), 
      shift1(0.0), shift2(0.0), specC(0), spec2rev(true)
  {
    // EMPTY
  }

  SpectrumPair::~SpectrumPair()
  {
    // EMPTY
  }

  SpectrumPair & SpectrumPair::operator=(const SpectrumPair & that)
  {
    spec1 = that.spec1;
    spec2 = that.spec2;
    shift1 = that.shift1;
    shift2 = that.shift2;
    score1 = that.score1;
    score2 = that.score2;
    spec2rev = that.spec2rev;
    return *this;
  }

  void SpectrumPair::output(ostream &out, char separator)
  {
    out << spec1 + 1 << separator << spec2 + 1 << separator << specC + 1
        << separator << shift1 << separator << shift2 << separator << score1
        << separator << score2 << endl;
  }

  void SpectrumPair::serialize(vector<float> &v)
  {
    v.resize(6);
    v[0] = spec1;
    v[1] = spec2;
    v[2] = shift1;
    v[3] = shift2;
    v[4] = score1;
    v[5] = score2;
  }

  unsigned int SpectrumPair::loadSz()
  {
    return 6;
  }

  void SpectrumPair::load(vector<float> &v)
  {
    if (v.size() < 6)
      return;
    spec1 = (unsigned int)v[0];
    spec2 = (unsigned int)v[1];
    shift1 = v[2];
    shift2 = v[3];
    score1 = v[4];
    score2 = v[5];
  }

} // namespace specnets
