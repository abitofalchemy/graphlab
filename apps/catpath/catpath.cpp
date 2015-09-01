/*-------*---------*---------*---------*---------*---------*---------*---------*
 * Author: s. aguinaga saguinag (at) nd dot edu
 *
 * catpath
 *
 * Star Log:
 * ---------
 * 11Aug15: Starting from scratch
 *
 * * This work is derive from:
 * Copyright (c) 2009 Carnegie Mellon University.
 *	  http://www.graphlab.ml.cmu.edu
 *
 */

#include <vector>
#include <string>
#include <fstream>
#include <limits.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/unordered_set.hpp>

#include <graphlab.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/macros_def.hpp>
#include <graphlab/ui/metrics_server.hpp>

#define COUT std::cout
#define ENDL std::endl


bool DIRECTED_SSSP = true;
bool PER_VERTEX_COUNT = false;

/**
 *  * \brief The type used to measure distances in the graph.
 *   */
typedef float distance_type ;
typedef   int wiki_page_type; // namespace


/**
 * * \ a custom struc we call Quad
 * * /
 * */
struct Quad {
	graphlab::vertex_id_type dest_art;
	distance_type from_last_art;
	distance_type from_src;
  graphlab::vertex_id_type last_cat;
  graphlab::vertex_id_type last_art;
	Quad(): dest_art(0) {}

	Quad(graphlab::vertex_id_type a, distance_type b, distance_type c, graphlab::vertex_id_type _last_art, graphlab::vertex_id_type d){
		dest_art	  = a;
		from_last_art = c;
		from_src	  = b;
	last_art = _last_art;
	last_cat = d;
	}

	void save(graphlab::oarchive& oarc) const {
	oarc << dest_art << from_last_art << from_src;}
	void load(graphlab::iarchive& iarc) {
	iarc >> dest_art >> from_last_art >> from_src;}

};

struct Triple: graphlab::IS_POD_TYPE
{
  distance_type from_last_art;
  distance_type from_src;
  graphlab::vertex_id_type last_cat;
  graphlab::vertex_id_type last_art;

  Triple():from_src(0){}
  Triple(distance_type b, distance_type c, graphlab::vertex_id_type _last_art, graphlab::vertex_id_type d)
  {
	from_last_art = c;
	from_src = b;
	last_cat = d;
	last_art = _last_art;
 	}
};


/**
 *  * \brief Get the other vertex in the edge.
 *   */
/**
 * \brief This class is used as the gather type.
 */
struct min_distance_type : graphlab::IS_POD_TYPE {
  distance_type dist;
  min_distance_type(distance_type dist =
				  std::numeric_limits<distance_type>::max()) : dist(dist) { }
  min_distance_type& operator+=(const min_distance_type& other) {
	dist = std::min(dist, other.dist);
	return *this;
	} // ends min_distance_type
};




/**
 *  * \brief The current distance of the vertex.
 *   */
struct vertex_data {
  distance_type  dist;
  wiki_page_type type; // 0 or 14 (namespace)
  boost::unordered_set<graphlab::vertex_id_type> vid_set;
  boost::unordered_set<graphlab::vertex_id_type> seen;
  graphlab::vertex_id_type prev_art; // +
  distance_type cat_dist_from_prev;
  std::list<Quad> msg_q;
  bool isDead;
  bool sent;
  vertex_data(wiki_page_type type = 0,
			  distance_type dist = std::numeric_limits<distance_type>::max(),
			  bool isDead=false,
			  bool sent=false) :
  dist(dist), type(type),isDead(isDead),sent(sent) { 		}


  void save(graphlab::oarchive& oarc) const {
	oarc << dist << type << isDead << sent;
	oarc << vid_set << seen ;
	uint32_t size = msg_q.size();
	oarc << size;
	for(std::list<Quad>::const_iterator it=msg_q.begin(); it != msg_q.end(); ++it ) {
	  oarc << *it;
	}
  }

  void load(graphlab::iarchive& iarc) {
	uint32_t msg_q_size;
	Quad tmp_quad;
	iarc >> dist >> type >> isDead >> sent;
	iarc >> vid_set >> seen;
	iarc >> msg_q_size;
	while(msg_q_size > 0) {
	  iarc >> tmp_quad;
	  msg_q.push_back(tmp_quad);
	  msg_q_size--;
	}
  }

};


/**
 *  * \brief The distance associated with the edge.
 *   */
struct edge_data{
  distance_type dist;
  edge_data() : dist(1) {}
  edge_data(distance_type _dist): dist(_dist) {}
  edge_data(min_distance_type _dist) {
	dist = _dist.dist;
  }
  void save(graphlab::oarchive& oarc) const { oarc << dist;}
  void load(graphlab::iarchive& iarc) { iarc >> dist;}
}; // end of edge data

struct our_msg {

  std::vector<Quad> collection;
  our_msg& operator+=(const our_msg& other){

  std::map<graphlab::vertex_id_type, Triple> temp;

  //COUT  << "om: other coll size: " << other.collection.size() << ENDL;
  //COUT  << "om: coll size: " << collection.size() << ENDL;
  for(unsigned int i=0; i< other.collection.size(); i++){
	Triple tp3(other.collection[i].from_src,
			   other.collection[i].from_last_art,
			   other.collection[i].last_art,
			   other.collection[i].last_cat);
	//std::cout << "om: Tp3:fs "<< other.collection[i].from_src	  << std::endl;
	//std::cout << "om: Tp3:fla "<< other.collection[i].from_last_art << std::endl;
	//COUT << "om: other coll dest_art: " <<other.collection[i].dest_art << ENDL;
	temp[other.collection[i].dest_art] = tp3;
  }

  for(unsigned int i=0; i< collection.size(); i++){
	if(temp.find(collection[i].dest_art) != temp.end()) {
	  if(temp[collection[i].dest_art].from_src < collection[i].from_src){

		collection[i].from_src = temp[collection[i].dest_art].from_src;
		collection[i].from_last_art = temp[collection[i].dest_art].from_last_art;
		collection[i].last_cat = temp[collection[i].dest_art].last_cat;
		collection[i].last_art = temp[collection[i].dest_art].last_art;
		temp[collection[i].dest_art].from_src = -2; // just to mark it
	  }
	}
  }

  std::map<graphlab::vertex_id_type, Triple>::iterator it;

  //put the ones not included in the collection
  for( it = temp.begin(); it!= temp.end(); it++){
	if( (*it).second.from_src != -2 ){
	  Quad tp2(it->first, it->second.from_src, it->second.from_last_art, it->second.last_art, it->second.last_cat);
	  collection.push_back(tp2);
	}
  }

	/*
	 std::set<graphlab::vertex_id_type> list_of_dest;

	 for(unsigned int i=0; i < collection.size(); i++)
	 {
	 list_of_dest.insert( collection[i].dest_art );
	 }


	 //Look for found articles that are searching for a new destination. If we
	 //find a message in search of an article that has already been found then
	 //we should kill that message.
	 for(std::vector<Quad>::iterator it=collection.begin(); it != collection.end();)
	 {
	 if( list_of_dest.find((*it).last_art) != list_of_dest.end()){
	 it = collection.erase(it);
	 }else{
	 ++it;
	 }
	 }

	 */

	return *this;
  }

  // serialize
  void save(graphlab::oarchive& oarc) const {
	oarc << collection;
	}

  // deserialize
  void load(graphlab::iarchive& iarc) {
	iarc >> collection;
		}
}; // ends our_msg



/**
 *  * \brief The graph type encodes the distances between vertices and
 *   * edges
 *	*/
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

inline graph_type::vertex_type
  get_other_vertex(const graph_type::edge_type& edge,
				   const graph_type::vertex_type& vertex) {
	return vertex.id() == edge.source().id()? edge.target() : edge.source();
}

struct set_union_gather {
  boost::unordered_set<graphlab::vertex_id_type> vid_set;
  /** Combining with another collection of vertices.
   * * Union it into the current set.
   ***/
  set_union_gather& operator+=(const set_union_gather& other) {
	foreach(graphlab::vertex_id_type othervid, other.vid_set) {
	  vid_set.insert(othervid);
	}
	return *this;
  }

  // serialize
  void save(graphlab::oarchive& oarc) const {
	oarc << vid_set;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
	iarc >> vid_set;
  }
};


class add_neighbours :
  public graphlab::ivertex_program<graph_type, set_union_gather>,
  /* I have no data. Just force it to POD */
  public graphlab::IS_POD_TYPE
{
  public:
  
  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
							 const vertex_type& vertex) const {
	return graphlab::IN_EDGES; // 08Aug15/SA: changed this to OUT
  }

  /*
   * For each edge, figure out the ID of the "other" vertex
   * and accumulate a set of the neighborhood vertex IDs.
   */
  gather_type gather(icontext_type& context,
					 const vertex_type& vertex,
					 edge_type& edge) const
  {
		set_union_gather gather;
		//std::cout<<vertex.id()<<" gather\n";
		// Insert the opposite end of the edge IF the opposite end has
		// ID greater than the current vertex
		// If we are getting per vertex counts, we need the entire neighborhood
		// if edge.source().id() == vertex.id()
		//	otherid = edge.target().id()
		// else
		//	otherid = edge.source().id()
    
    //if (vertex.id() == 13)
    //  COUT << " "<< vertex.id()  << ENDL;
    vertex_id_type otherid = edge.source().id();
    //edge.source().id() == vertex.id() ?
                             //edge.target().id() : edge.source().id();
    
    
		if (vertex.data().type== 0 && (edge.source().data().type == 0 ||
                                   edge.source().data().type == 14) && vertex.num_in_edges()>1)
		{
		  gather.vid_set.insert(otherid);
		  // std::cout << "  "<< vertex.id() << ", nei (other): "<< otherid << std::endl;
      //COUT << " "<< edge.source().id() << " "<< edge.target().id()  << ENDL;
		}
    
    return gather;
  }

  /*
   * the gather result now contains the vertex IDs in the neighborhood.
   * store it on the vertex.
   */
  void apply(	icontext_type& context,
			 vertex_type& vertex,
			 const gather_type& neighborhood) {
		vertex.data().vid_set = neighborhood.vid_set;
  } // end of apply

  /*
   * Scatter over all edges to compute the intersection.
   * I only need to touch each edge once, so if I scatter just on the
   * out edges, that is sufficient.
   */
  edge_dir_type scatter_edges(icontext_type& context,
							  const vertex_type& vertex) const {
	//std::cout << graphlab::NO_EDGES << std::endl;
	return graphlab::NO_EDGES; // SA: changed this to OUT
  }

  // void scatter(icontext_type& context,
  //				 const vertex_type& vertex,
  //				 edge_type& edge) const {
  //   not needed
  //	 }




};// ends add_neighbours


/**
 * \brief The single source shortest path vertex program.
 */
class main_algo : public graphlab::ivertex_program<graph_type,
				  graphlab::empty,
				  our_msg>
{
  std::vector<Quad> vrtx_dists;
  public:
	void init(icontext_type& context,
			  const vertex_type& vertex,
			  const our_msg& msg) { vrtx_dists = msg.collection; }

  /**
   * \ brief We use the messaging model to compute the SSSP update
   * / Gather on NO edges?
   */
  edge_dir_type gather_edges(icontext_type& context,
							 const vertex_type& vertex) const {
	return graphlab::NO_EDGES;
  }// end of gather_edges

  /**
   * \brief The apply phase: If the distance is smaller then update
   */
  void apply(icontext_type&		 context,
						 vertex_type&		   vertex,
			 const graphlab::empty& empty)
  {
  	//COUT << "apply phase: " << vertex.id() << ENDL;
		distance_type tp = std::numeric_limits<distance_type>::max(); // inf
		distance_type cat_dist_from_prev = 0;
		graphlab::vertex_id_type prev_art = 0;
    
		if(vertex.data().sent == true){
	  	vertex.data().isDead = true;
		} else if(vertex.data().type==0){ // if a article **************************
		  for(unsigned int i=0; i<vrtx_dists.size();i++){
				// From our_msg
				//gets minimum distance from incoming category edges
				if(vrtx_dists[i].from_src < tp){
				  tp = vrtx_dists[i].from_src;
				  prev_art = vrtx_dists[i].last_art;
				  cat_dist_from_prev = vrtx_dists[i].from_last_art;
				} // ends if
			} // ends for loop

		  // Is the minimum we just found, less than the distance we currently
		  // store? If so, store new minimum.
		  if( tp < vertex.data().dist ){
				vertex.data().dist	 = tp;
				vertex.data().prev_art = prev_art;
				vertex.data().cat_dist_from_prev = cat_dist_from_prev;
				vertex.data().sent = true; // mark the node as seen
			}
			//COUT << "apply: v.dist= " << vertex.data().dist << ENDL;
		} //end of if an article

		else if(vertex.data().type==14){ // if a category -- very unclear ********
		  std::map<graphlab::vertex_id_type, Triple> temp;

		  for(unsigned int i=0; i < vrtx_dists.size(); i++)	{
				//COUT << "  cat|vrtx_dists[i].dest_art "<<vrtx_dists[i].dest_art<<ENDL;
				vertex.data().seen.insert(vrtx_dists[i].dest_art);
				vrtx_dists[i].from_last_art +=1;
				vrtx_dists[i].from_src +=1;
				if(vertex.data().msg_q.size() > 0){
				  vertex.data().msg_q.pop_front(); //removed in the previous scatter
				}

				vertex.data().msg_q.push_back( vrtx_dists[i] ) ;
			}// ends for

		} // ends if|else


  } // completes apply

  /**
   * \brief Determine if SSSP should run on all edges or just in edges
   */
  edge_dir_type scatter_edges(icontext_type& context,
							  const vertex_type& vertex) const {
	return graphlab::OUT_EDGES; // we scatter over out-going edges
  } // end of scatter_edges

  /**
   * \brief The scatter function just signal adjacent pages
   */
  void scatter( icontext_type& context,
								const vertex_type& vertex,
			   				edge_type& edge) const
  {
		const vertex_type other = get_other_vertex(edge, vertex);
		//COUT<< "scatter: " << vertex.id() << "->" << other.id() << ENDL;

		if(vertex.data().type == 14 && other.data().type==14) //category to category
		{
		  our_msg for_cat;
		  Quad msg = vertex.data().msg_q.front(); // get the message
		  if(other.data().seen.find(msg.dest_art) /* if cat node has seen */
				 != other.data().seen.end())
		  {
				msg.from_src += 1;
				msg.from_last_art += 1;
				msg.last_cat = vertex.id();
				for_cat.collection.push_back(msg); //updating distance
		  }

	  if(for_cat.collection.size() !=0){
		context.signal(other, for_cat);
	  }

	} // ends cat to cat
	else if (	vertex.data().type  == 14 &&
			 				other.data().type   == 0 &&
			 				other.data().isDead == false) //category to article ***************
	{
	  our_msg for_art;
	  for(unsigned int i=0; i<vrtx_dists.size(); i++){
			// If it exists in the add_neighbours: print out the dest art
//			COUT 	<< " scatter: (v,o) " << vertex.id() << "," << other.id()
//				 		<< "   dest art: " << vrtx_dists[i].dest_art << std::endl;

			if(other.data().vid_set.find( vrtx_dists[i].dest_art ) !=
			   other.data().vid_set.end() )
			{
//			  printf("   cat->nei: %llu->%llu\n", vrtx_dists[i].dest_art, other.id());

			  Quad tp2(vrtx_dists[i].dest_art,
					   vrtx_dists[i].from_src+1,
					   vrtx_dists[i].from_last_art+1,
					   vrtx_dists[i].last_art,
					   vertex.id());
			  if(vrtx_dists[i].from_src != std::numeric_limits<distance_type>::max())
			  {
				for_art.collection.push_back(tp2);
			  }
			}//ends if
	  } // ends for
	  if(for_art.collection.size() !=0){
			// signals neighbors ?
			context.signal(other,for_art);
	  	}
	} // ends else if
	else if(vertex.data().type == 0 &&
					other.data().type  == 14 &&
				  vertex.data().isDead == false)//article to category
	{
	  //COUT << " art->cat: " << vertex.id() << " -> " << other.id() << ENDL;
	  our_msg for_cat1;
	  Quad tp2(other.id(), vertex.data().dist, 1, vertex.id(), 0); // changed from other v
	  for_cat1.collection.push_back(tp2);

	  if(vertex.data().dist != std::numeric_limits<distance_type>::max())
		context.signal(other, for_cat1);
	}

  } // end of scatter

  // serialize
  void save(graphlab::oarchive& oarc) const {	oarc << vrtx_dists; }

  // deserialize
  void load(graphlab::iarchive& iarc) { iarc >> vrtx_dists;  }

}; // end of shortest path vertex program


/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
//------Output of the final graph----------
struct shortest_path_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
	std::stringstream strm;

	if(vtx.data().dist != std::numeric_limits<distance_type>::max() )
	  strm << vtx.id() << "\t" << "ns:"
		   << "\t" << vtx.data().type
		   << "\t" <<vtx.data().dist
		   << "\t" <<vtx.data().prev_art
		   << "\t" <<vtx.data().cat_dist_from_prev<<"\n";

	//	 boost::unordered_set<graphlab::vertex_id_type>::iterator it;
	// 	for( it = vtx.data().vid_set.begin(); it !=vtx.data().vid_set.end(); it++)
	//		 strm<<*it<<" ";


	return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of shortest_path_writer


bool line_parser_art(graph_type& graph,
					 const std::string& filename,
					 const std::string& textline)
{
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
  graphlab::vertex_id_type tid;
  // 	wiki_page_type type =1;

  // first entry in the line is a vertex ID
  strm >> vid;
  strm >> tid;

  // //std::cout<<textline<<std::endl;
  //dc.cout()<<textline<<std::endl;
  if(tid!=vid){

	//	//std::cout<<"one \n"
	//	//std::cout<<"adding edges \n";
	graph.add_edge(vid, tid);
	////std::cout<<"two \n";
  }
  return true;
}

bool line_parser_categ (graph_type& graph,
						const std::string& filename,
						const std::string& textline)
{
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
  graphlab::vertex_id_type tid;
  //   wiki_page_type type =1;
  // first entry in the line is a vertex ID
  strm >> vid;
  strm >> tid;

  // //std::cout<<textline<<std::endl;
  //dc.cout()<<textline<<std::endl;
  if(tid!=vid){
	//	//std::cout<<"one \n"
	//	//std::cout<<"adding edges \n";
	graph.add_edge(vid, tid);
	graph.add_edge(tid, vid);
	////std::cout<<"two \n";
  }
  return true;
}
bool all_vertex_parser(graph_type& graph,
					   const std::string& filename,
					   const std::string& textline)
{
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
  wiki_page_type type;

  // first entry in the line is a vertex ID
  strm >> vid;
  strm >> type;

  // //std::cout<<textline<<std::endl;
  if(type==0 || type ==14)
	graph.add_vertex(vid, vertex_data(type));
  return true;

}

/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
struct graph_writer {
//  std::string save_vertex(const graph_type::vertex_type& vtx) {
//	std::stringstream strm;
//	strm << vtx.id() << "\t" << vtx.data().dist << "\n";
//	return strm.str();
//  }
//  std::string save_edge(graph_type::edge_type e) { return ""; }
  std::string save_vertex(const graph_type::vertex_type& vtx) {
	return "";
	}
  std::string save_edge(graph_type::edge_type e) {

	std::stringstream strm;
	strm<<e.source().data().type<<" "<<e.target().data().type<<"\n";
	return strm.str();
  }
}; // end of graph_writer


int main( int argc, char** argv) {
  // Initialization block
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options
  graphlab::command_line_options clopts("Catpath: A Weninger SP Algorithm.");
  std::string graph_dir;
  std::string format = "snap";
  std::string exec_type = "synchronous";
  size_t powerlaw = 0;

  //std::vector<graphlab::vertex_id_type> sources;
  std::vector<unsigned int> sources;
  bool max_degree_source = false;
  clopts.attach_option( "graph", graph_dir,
												"The graph file.  If none is provided "
												"then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option( "format", format,  "graph format");
  clopts.attach_option( "source", sources, "The source vertices");
  clopts.attach_option( "max_degree_source", max_degree_source,
												"Add the vertex with maximum degree as a source");
  clopts.add_positional("source");
  clopts.attach_option( "directed", DIRECTED_SSSP, "Treat edges as directed.");
  clopts.attach_option( "engine", exec_type,
												"The engine type synchronous or asynchronous");
  clopts.attach_option( "powerlaw", powerlaw,
												"Generate a synthetic powerlaw out-degree graph. ");
  std::string saveprefix;
  clopts.attach_option( "saveprefix", saveprefix,
												"If set, will save the resultant catpath to a "
												"sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
		dc.cout() << "Error in parsing command line arguments." << std::endl;
		return EXIT_FAILURE;
  }

  // Begin
  // Build the graph
  graph_type graph(dc, clopts);
  graph.load("/Users/saguinag/Research/datasets/enwiki/page.txt" ,
			   all_vertex_parser);
  graph.load("/Users/saguinag/Research/datasets/page_links_topranked.tsv",
			   line_parser_art);
  graph.load("/Users/saguinag/Research/datasets/enwiki/catlink.txt",
			   line_parser_categ);
  logstream(LOG_INFO) << "finish loading all graphs\n";


  // Must call finalize before querying the graph
  graph.finalize();
  //    dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
  //              << "#edges:	 " << graph.num_edges() << std::endl;

  // Source nodes
  if(sources.empty()) {
	if (max_degree_source == false) {
		dc.cout()<< "No source vertex provided. Add an argument specifying the"
					"source id (e.g., --source=# )"
		<< std::endl;
		sources.clear();
		sources.push_back(303);//abort();
	  }
  }

  // Running The Engines -------------------------------------------------------
  // Class `add_neighbours` collects each nodes' out-going neighbors into a set
  graphlab::omni_engine<add_neighbours> engine(dc, graph, exec_type, clopts);
  engine.signal_all();
  engine.start();
  
  //abort();

  // Class main_algo computes the catpaths
  graphlab::omni_engine<main_algo> engine2(dc, graph, exec_type, clopts);
  our_msg init_msg;
  Quad seed_values(sources[0], /*dest_art*/
				   0, /*dist from_last_art*/
				   0, /*dist from_src*/
				   0, /*v.id last_cat*/
				   0  /*v.id last_art*/);
  init_msg.collection.push_back(seed_values);
  engine2.signal(sources[0], init_msg); // starting node and the msg.
  engine2.start();

  const float runtime = engine2.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime
						<< " seconds." << std::endl;

  // Save the final graph -----------------------------------------------------
  //isaveprefix = "/data/saguinag/Results/catpath_";
  saveprefix = "/tmp/catpath_";
  if (saveprefix != "") {
	graph.save(saveprefix, shortest_path_writer(),
			   false,	// do not gzip
			   true,	 // save vertices
			   false);   // do not save edges
  }

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
} // end of main
