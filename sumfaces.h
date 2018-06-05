#include "freespace.h"

class Feature : public RefCnt {
 public:
  Points p;
  Feature(PTR<Point> a) { p.push_back(a); }
  Feature(PTR<Point> a, PTR<Point> b) { p.push_back(a); p.push_back(b); }
  Feature(PTR<Point> a, PTR<Point> b, PTR<Point> c) { p.push_back(a); p.push_back(b); p.push_back(c); }
  Feature(std::vector<PTR<Point> >& p) { this->p.insert(this->p.end(), p.begin(), p.end()); }

  PTR<Feature> rotate(PTR<Point> sin_cos_alpha) {
    Feature * f = new Feature(p);
    for (int i=0; i<f->p.size(); i++) {
      f->p[i] = new RotationPoint(f->p[i], sin_cos_alpha);
    }
    return f;
  }

  PTR<Feature> sum(PTR<Feature> that) {
    std:std::vector<PTR<Point> > p_new;
    if (that->p.size() == 1) {
      for(int i=0; i<this->p.size(); i++)
        p_new.push_back(new SumPoint(that->p[0], this->p[i]));
      return new Feature(p_new);
    }

    if (this->p.size() == 1) {
      for(int i=0; i<that->p.size(); i++)
        p_new.push_back(new SumPoint(this->p[0], that->p[i]));
      return new Feature(p_new);
    }

    if (this->p.size() == 2 && this->p.size() == 2) {
      p_new.push_back(new SumPoint(this->p[0], that->p[0]));
      p_new.push_back(new SumPoint(this->p[0], that->p[1]));
      p_new.push_back(new SumPoint(this->p[1], that->p[1]));
      p_new.push_back(new SumPoint(this->p[1], that->p[0]));
      return new Feature(p_new);
    }

    assert(false);
  }
};

bool * candidatePairs(PTR<Feature> v, PTR<Feature> f, PTR<Point> sin_cos_alpha, PTR<Point> sin_cos_alpha_sample, std::vector<Polyhedron*> & blockspaces, int n_samples);
