//Adapted from Lindsey Gray
//https://github.com/lgray/cmssw/blob/new_hgc_clusters_inprogress/RecoParticleFlow/PFClusterProducer/src/HGCClusterizer.cc#L369-L539

#ifndef Clusterizer_hh
#define Clusterizer_hh

#include "HGCSSRecoHit.hh"
#include "HGCSSCluster.hh"
// helpful tools
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include <unordered_set>

//local tools
namespace {
  std::pair<float,float> minmax(const float a, const float b) {
    return ( b < a ? 
	     std::pair<float,float>(b, a) : 
	     std::pair<float,float>(a, b)   );
  }

  template<typename T>
  KDTreeCube fill_and_bound_kd_tree(const std::vector<T>& points,
				    const std::vector<bool>& usable,
				    std::vector<KDTreeNodeInfoT<unsigned,3> >& nodes) {
    std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
    for( unsigned i = 0 ; i < points.size(); ++i ) {
      if( !usable[i] ) continue;
      const auto& pos = points[i].position();
      nodes.emplace_back(i, (float)pos.X(), (float)pos.Y(), (float)pos.Z());

      //std::cout << " ele " << i << " position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << " nodes size:" << nodes.size() << std::endl;

      if( i == 0 ) {
	minpos[0] = pos.X(); minpos[1] = pos.Y(); minpos[2] = pos.Z();
	maxpos[0] = pos.X(); maxpos[1] = pos.Y(); maxpos[2] = pos.Z();
      } else {
	minpos[0] = std::min((float)pos.X(),minpos[0]);
	minpos[1] = std::min((float)pos.Y(),minpos[1]);
	minpos[2] = std::min((float)pos.Z(),minpos[2]);
	maxpos[0] = std::max((float)pos.X(),maxpos[0]);
	maxpos[1] = std::max((float)pos.Y(),maxpos[1]);
	maxpos[2] = std::max((float)pos.Z(),maxpos[2]);
      }
    }

    //std::cout << " min and max X = " << minpos[0] << " " << maxpos[0] << std::endl
    //<< " min and max Y = " << minpos[1] << " " << maxpos[1] << std::endl
    //<< " min and max Z = " << minpos[2] << " " << maxpos[2] << std::endl;

    return KDTreeCube(minpos[0],maxpos[0],
		      minpos[1],maxpos[1],
		      minpos[2],maxpos[2]);
  }

  bool greaterByEnergy(const std::pair<unsigned,double>& a,
		       const std::pair<unsigned,double>& b) {
    return a.second > b.second;
  }

  class QuickUnion{    
  public:
    QuickUnion(const unsigned NBranches) {
      _count = NBranches;
      _id.resize(NBranches);
      _size.resize(NBranches);
      for( unsigned i = 0; i < NBranches; ++i ) {
	_id[i] = i;
	_size[i] = 1;
      }
    }
    
    int count() const { return _count; }
    
    unsigned find(unsigned p) {
      while( p != _id[p] ) {
	_id[p] = _id[_id[p]];
	p = _id[p];
      }
      return p;
    }
    
    bool connected(unsigned p, unsigned q) { return find(p) == find(q); }
    
    void unite(unsigned p, unsigned q) {
      unsigned rootP = find(p);
      unsigned rootQ = find(q);
      _id[p] = q;
      
      if(_size[rootP] < _size[rootQ] ) { 
	_id[rootP] = rootQ; _size[rootQ] += _size[rootP]; 
      } else { 
	_id[rootQ] = rootP; _size[rootP] += _size[rootQ]; 
      }
      --_count;
    }
    std::vector<unsigned> _id;
    std::vector<unsigned> _size;
    int _count;
    
  };//class
}//namespace

class Clusterizer{
  typedef KDTreeLinkerAlgo<unsigned,3> KDTree;
  typedef KDTreeNodeInfoT<unsigned,3> KDNode;
  typedef std::pair<unsigned,unsigned> HitLink;
  typedef std::unordered_multimap<unsigned,unsigned> LinkMap;
  typedef std::unordered_set<unsigned> UniqueIndices;

public:
  Clusterizer(unsigned debug=0);
  ~Clusterizer();  
  void buildClusters(std::vector<HGCSSRecoHit> *rechitvec,
		     const std::vector<bool>&,
		     const std::vector<bool>&, 
		     HGCSSClusterVec &);
 
private:

  void build2DCluster(const std::vector<HGCSSRecoHit>  & rechitvec,
		      const std::vector<bool>& rechitMask,
		      const std::vector<bool>& seedable,
		      const unsigned current_index,
		      std::vector<bool>& usable,
		      HGCSSCluster & cluster);

  void linkClustersInLayer(const HGCSSClusterVec & input_clusters,
			   HGCSSClusterVec & output);


  double _moliR;
  unsigned debug_;

  // used for rechit searching 
  std::vector<KDNode> _cluster_nodes, _hit_nodes, _found;
  KDTree _cluster_kdtree,_hit_kdtree;


};//end class


#endif
