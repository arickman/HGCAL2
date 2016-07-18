#include <iomanip>
#include <cmath>
// helpful tools
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include <unordered_set>

#include "Clusterizer.hh"

#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

Clusterizer::Clusterizer(unsigned debug){

  _moliR = 2.7;//cm
  debug_ = debug;

  // clean initial state for searching
  _cluster_nodes.clear(); _found.clear(); _cluster_kdtree.clear();
  _hit_nodes.clear(); _hit_kdtree.clear();

}

Clusterizer::~Clusterizer(){
}

void Clusterizer::
buildClusters(std::vector<HGCSSRecoHit> *rechitvec,
	      const std::vector<bool>& rechitMask,
	      const std::vector<bool>& seedable,
	      HGCSSClusterVec & output) {
  
  const std::vector<HGCSSRecoHit> & rechits = *rechitvec;
  std::vector<bool> usable_rechits(rechits.size(),true);
  std::unordered_set<double> unique_depths;

  HGCSSClusterVec clusters_per_layer;
  
  if (debug_) std::cout << " -- Building clusters out of " << rechits.size() << " hits: " << std::endl;

  // sort seeds by energy
  std::vector<unsigned> seeds;
  for( unsigned i = 0; i < rechits.size(); ++i ) {
    if( seedable[i] ) {
      auto pos = std::lower_bound(seeds.begin(),seeds.end(),i,
				  [&](const unsigned i, const unsigned j) {
				    return ( rechits[i].energy() >= 
					     rechits[j].energy()   );
				  });
      seeds.insert(pos,i);
    }
  }
  
  if (debug_) std::cout << " -- Size of seed vec: " << seeds.size() << std::endl;

  if (debug_>1){ 
    for ( unsigned i = 0; i < seeds.size(); ++i ) {
      std::cout << " Seed " << i << " index " << seeds[i] << " energy " << rechits[seeds[i]].energy() << std::endl;
    }
  }

  // get ready for initial topo clustering
  KDTreeCube kd_boundingregion = 
    fill_and_bound_kd_tree(rechits,usable_rechits,_hit_nodes);
  _hit_kdtree.build(_hit_nodes,kd_boundingregion);
  _hit_nodes.clear();
  // make topo-clusters that require that the energy goes
  // down with respect to the last rechit encountered
  // rechits clustered this way are locked from further use  
  for( const unsigned i : seeds ) {
    const auto& hit = rechits[i];
    if (debug_>1) std::cout << " - Seed idx " << i << " energy " << hit.energy() << " layer " << hit.layer() << std::endl;
    //if(hit.neighbours8().size() > 0 ) {
    HGCSSCluster layer_cluster;
    build2DCluster(rechits, rechitMask, seedable,
		   i, usable_rechits, 
		   layer_cluster);

    if (debug_>1) std::cout << " --- Cluster has : " << layer_cluster.recHitFractions().size() << " rechits associated. Energy of first ele: " << layer_cluster.recHitFractions().begin()->first->energy() << " weight = " << layer_cluster.recHitFractions().begin()->second << std::endl;

    unique_depths.insert(std::abs(hit.position().Z()));

    const auto& hAndFs = layer_cluster.recHitFractions();
    if( hAndFs.size() > 1 ) {
      layer_cluster.setLayer(hit.layer());
      layer_cluster.setSeed(hit.position());
      layer_cluster.setSeedEnergy(hit.energy());
      layer_cluster.calculatePosition();
      clusters_per_layer.push_back(std::move(layer_cluster));
    } else if ( hAndFs.size() == 1 ) {
      usable_rechits[i] = true;
    }
    
  }

  if (debug_) std::cout << " -- Number of clusters per layer: " << clusters_per_layer.size() << std::endl;

  _hit_kdtree.clear();
  
  HGCSSClusterVec z_linked_clusters;
  // use topo clusters to link in z
  linkClustersInLayer(clusters_per_layer,z_linked_clusters); 

  if (debug_) std::cout << " -- Number of clusters after linking in z: " << z_linked_clusters.size() << std::endl;

  std::vector<bool> usable_clusters(z_linked_clusters.size(),true);
  
  // stuff usable clusters into the output list
  for( unsigned i = 0; i < z_linked_clusters.size(); ++i ) {
    if( i >= usable_clusters.size() ) {
      output.push_back(z_linked_clusters[i]);
    } else if( usable_clusters[i] ) {
      output.push_back(z_linked_clusters[i]);
    }
  }
}//buildCluster


void Clusterizer::
build2DCluster(const std::vector<HGCSSRecoHit>  & rechitvec,
	       const std::vector<bool>& rechitMask,
	       const std::vector<bool>& seedable,
	       const unsigned current_index,
	       std::vector<bool>& usable,
	       HGCSSCluster & cluster){

  usable[current_index] = false;
  const HGCSSRecoHit & current_cell = rechitvec[current_index];

  if (debug_>1) std::cout << " -- Current index = " << current_index << " hit energy = " << current_cell.energy() << std::endl;

  cluster.addRecHitFraction(std::pair<HGCSSRecoHit*,double>(const_cast<HGCSSRecoHit*>(&current_cell),1.0));

  if (debug_>2) std::cout << " Cluster now has : " << cluster.recHitFractions().size() << " rechits associated. Energy of first ele: " << cluster.recHitFractions().begin()->first->energy() << std::endl;
  

  //CAMM: mm??
  double moliere_radius = _moliR;
  const ROOT::Math::XYZPoint pos = current_cell.position();
  
  auto x_rh = minmax(pos.x()+moliere_radius,pos.x()-moliere_radius);
  auto y_rh = minmax(pos.y()+moliere_radius,pos.y()-moliere_radius);
  //CAMM 1um ?? Need to change to 300um :/
  auto z_rh = minmax(pos.z()+0.03,pos.z()-0.03);

  KDTreeCube hit_searchcube((float)x_rh.first,(float)x_rh.second,
			    (float)y_rh.first,(float)y_rh.second,
			    (float)z_rh.first,(float)z_rh.second);
  std::vector<KDNode> found;
  _hit_kdtree.search(hit_searchcube,found);

  if (debug_>1) std::cout << " -- Number of closest neighbours found: " << found.size() << std::endl; 

  for( const KDNode& nbourpoint :found ) {
    // only cluster if not a seed, not used, and energy less than present
    if (debug_>2) std::cout << " Neighbour index: " << nbourpoint.data << " usable = " << usable[nbourpoint.data] << " seedable " << seedable[nbourpoint.data] << " energy " << rechitvec[nbourpoint.data].energy() << " rechitmask " << rechitMask[nbourpoint.data] << std::endl;
    const HGCSSRecoHit& nbour = rechitvec[nbourpoint.data];
    if( usable[nbourpoint.data] && !seedable[nbourpoint.data] &&
	nbour.energy() <= current_cell.energy() && // <= takes care of MIP sea
	rechitMask[nbourpoint.data]) {
      //std::cout << " search for next neighbours!" << std::endl;
      build2DCluster(rechitvec,rechitMask,seedable,nbourpoint.data,usable,cluster);
    }
    //else std::cout << " -- going to next found node." << std::endl;
  }
}//build2Dcluster



void Clusterizer::
linkClustersInLayer(const HGCSSClusterVec & input_clusters,
		    HGCSSClusterVec & output) {

  std::vector<bool> dummy(input_clusters.size(),true);
  KDTreeCube kd_boundingregion =
    fill_and_bound_kd_tree(input_clusters,dummy,_cluster_nodes);
  _cluster_kdtree.build(_cluster_nodes,kd_boundingregion);
  _cluster_nodes.clear();
  //const float moliere_radius2 = std::pow(moliere_radius,2.0);
  LinkMap back_links;
  // now link all clusters with in moliere radius for EE + HEF
  for( unsigned i = 0; i < input_clusters.size(); ++i ) {
    float moliere_radius = _moliR;

    const auto& incluster = input_clusters[i];
    const auto& pos = incluster.position();
    auto x = minmax(pos.X()+moliere_radius,pos.X()-moliere_radius);
    auto y = minmax(pos.Y()+moliere_radius,pos.Y()-moliere_radius);
    //CAMM: why 2* ?
    auto z = minmax(pos.Z()+2*moliere_radius,pos.Z()-2*moliere_radius);
    KDTreeCube kd_searchcube((float)x.first,(float)x.second,
			     (float)y.first,(float)y.second,
			     (float)z.first,(float)z.second);
    _cluster_kdtree.search(kd_searchcube,_found);
    for( const auto& found_node : _found ) {
      const auto& found_clus = input_clusters[found_node.data];
      const auto& found_pos = found_clus.position();
      const auto& diff_pos = found_pos - pos;
      //CAMM check 0.001 val
      if( diff_pos.rho() < moliere_radius && std::abs(diff_pos.Z()) > 0.03 ) {
	if( pos.mag2() > found_pos.mag2() ) {
	  back_links.insert(std::pair<unsigned,unsigned>(i,found_node.data));
	} else {
	  back_links.insert(std::pair<unsigned,unsigned>(found_node.data,i));
	}
      }
    }
    _found.clear();
  }
  // using back-links , use simple metric for now to get something working
  QuickUnion qu(input_clusters.size());
  unsigned best_match;
  float min_parameter;
  for( unsigned i = 0; i < input_clusters.size(); ++i ) {
    const auto& pos = input_clusters[i].position();
    const auto clusrange = back_links.equal_range(i);
    min_parameter = std::numeric_limits<float>::max();
    best_match = std::numeric_limits<unsigned>::max();
    for( auto connected = clusrange.first;
	 connected != clusrange.second; ++connected ) {
      const auto& pos_connected = input_clusters[connected->second].position();
      float angle = (pos_connected - pos).theta();
      if( pos.z() < 0.0f ) angle += M_PI;
      while( angle > M_PI ) angle -= 2*M_PI;
      while( angle < -M_PI ) angle += 2*M_PI;
      angle = std::abs(angle) + 0.001f;
      const float dist2 = (pos_connected - pos).Mag2();
      const float parm = dist2*angle*angle;
      if( parm < min_parameter ) {
	best_match = connected->second;
	min_parameter = parm;
      }
    }
    if( best_match != std::numeric_limits<unsigned>::max() ) {
      qu.unite(i,best_match);
    }
  }
  LinkMap merged_clusters;
  UniqueIndices roots;
  for( unsigned i = 0; i < input_clusters.size(); ++i ) {
    const unsigned root = qu.find(i);
    roots.insert(root);
    merged_clusters.insert(std::pair<unsigned,unsigned>(root,i));
  }
  //std::cout << roots.size() << " final clusters!" << std::endl;

  
  unsigned iclus = 0;
  for( const auto& root : roots ) {
    HGCSSCluster merged_cluster;
    float max_energy = 0;
    HGCSSRecoHit* seed_hit = 0;
    auto range = merged_clusters.equal_range(root);
    for( auto clus = range.first; clus != range.second; ++clus ) {

      const auto& hAndFs = input_clusters[clus->second].recHitFractions();
      for( const auto& hAndF : hAndFs ) {
	merged_cluster.addRecHitFraction(hAndF);
	if( hAndF.first->energy() > max_energy ) {
	  max_energy = hAndF.first->energy();
	  seed_hit = hAndF.first;
	}
      }
    }
    if (!seed_hit){
      std::cout << " Problem! Seed hit not found." << std::endl;
      exit(1);
    }
    merged_cluster.setSeed(seed_hit->position());
    merged_cluster.setSeedEnergy(max_energy);
    merged_cluster.setLayer(seed_hit->layer());
    merged_cluster.calculatePosition();
    output.push_back(merged_cluster);
    ++iclus;
  }
  
  
  _found.clear();
  _cluster_kdtree.clear();
  
}//link clusters
