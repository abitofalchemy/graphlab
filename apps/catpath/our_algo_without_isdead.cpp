/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
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

/**
 * \brief The type used to measure distances in the graph.
 */
typedef float distance_type;

typedef int		wiki_page_type;

/**
 * \brief The current distance of the vertex.
 */
struct vertex_data {
  distance_type dist;
	wiki_page_type type;
	boost::unordered_set<graphlab::vertex_id_type> vid_set;

  vertex_data(wiki_page_type type = 0, distance_type dist = std::numeric_limits<distance_type>::max()) :
    dist(dist), type(type) { 		}

	void save(graphlab::oarchive& oarc) const { 
		oarc << dist << type << vid_set; 
	}
	
	void load(graphlab::iarchive& iarc) { 
		iarc >> dist >> type >> vid_set; 
	} 

};

/**
 * \brief The distance associated with the edge.
 */
struct edge_data : graphlab::IS_POD_TYPE {
  distance_type dist;
  edge_data(distance_type dist = 1) : dist(dist) { }
}; // end of edge data


/**
 * \brief The graph type encodes the distances between vertices and
 * edges
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief Get the other vertex in the edge.
 */
inline graph_type::vertex_type
get_other_vertex(const graph_type::edge_type& edge,
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}


/**
 * \brief Use directed or undireced edges.
 */
bool DIRECTED_SSSP = false;


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
  }
};

bool PER_VERTEX_COUNT = false;

struct triple : graphlab::IS_POD_TYPE {
	graphlab::vertex_id_type vid;
	distance_type from_last_art;
	distance_type from_src;

	triple(): vid(0) {}

	triple(graphlab::vertex_id_type a, distance_type b, distance_type c){
		vid =a;
		from_last_art = c;
		from_src = b;
	}

};

struct set_union_gather {
	  boost::unordered_set<graphlab::vertex_id_type> vid_set;

		  /*
			 * * Combining with another collection of vertices.
			 * * Union it into the current set.
			 * */
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
													                   
	
struct our_msg {

std::vector<triple> collection;

our_msg& operator+=(const our_msg& other){
	

	std::map<graphlab::vertex_id_type, std::pair<distance_type,distance_type> > temp;
	
	for(int i=0; i< other.collection.size(); i++){
   temp[other.collection[i].vid]=std::make_pair(other.collection[i].from_src,other.collection[i].from_last_art);
	}
	for(int i=0; i< collection.size(); i++){
		if(temp.find(collection[i].vid) != temp.end()){
     if(temp[collection[i].vid].first< collection[i].from_src){
        collection[i].from_src= temp[collection[i].vid].first;
		    collection[i].from_last_art = temp[collection[i].vid].second;
				temp[collection[i].vid].first = -2; // just to mark it
	  	}
		}
	}
	std::map<graphlab::vertex_id_type,std::pair<distance_type,distance_type> >::iterator it;

	//put the ones not icnluded in the collection 
	for( it = temp.begin(); it!= temp.end(); it++){
    if( (*it).second.first != -2 ){
			triple tp2(it->first,(*it).second.first,(*it).second.second);
		collection.push_back(tp2);
		}	
	}
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
};	

class add_neighbours :
      public graphlab::ivertex_program<graph_type,
                                      set_union_gather>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE {
public:
  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::IN_EDGES;
  }

  /*
* For each edge, figure out the ID of the "other" vertex
* and accumulate a set of the neighborhood vertex IDs.
*/
  gather_type gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    set_union_gather gather;
    //std::cout<<vertex.id()<<" gather\n";
    // Insert the opposite end of the edge IF the opposite end has
    // ID greater than the current vertex
    // If we are getting per vertex counts, we need the entire neighborhood
    vertex_id_type otherid = edge.source().id() == vertex.id() ?
                             edge.target().id() : edge.source().id();
     if(vertex.data().type== 0 && edge.target().data().type == 0)
		 gather.vid_set.insert(otherid);
    return gather;
  }

  /*
* the gather result now contains the vertex IDs in the neighborhood.
* store it on the vertex.
*/
  void apply(icontext_type& context, vertex_type& vertex,
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
    return graphlab::NO_EDGES;
  }
  
  // void scatter(icontext_type& context,
//                 const vertex_type& vertex,
//                 edge_type& edge) const {
//   not needed
//     }



};

class main_algo:
      public graphlab::ivertex_program<graph_type,
																			graphlab::empty, 
                                      our_msg>  // our_msg is my message type not gather type 
    
       {

 std::vector<triple> temp1;
public:
  void init(icontext_type& context, const vertex_type& vertex,
           const our_msg& msg) {
    temp1 = msg.collection;
  } 


  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  void apply(icontext_type& context, vertex_type& vertex,
              const graphlab::empty& empty) {
    distance_type tp=INT_MAX;
  if(vertex.data().type==0){ // if a article
		for(int i=0; i<temp1.size();i++){
     if(temp1[i].from_src < tp){
			 tp=temp1[i].from_src;
	  	}
    }
		if( tp< vertex.data().dist){
			vertex.data().dist = tp;
    }
	}//end of if an article	
	 else if(vertex.data().type==14){//its a category
    // dont do anything

	 }
  } // end of apply

  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }
  
  // void scatter(icontext_type& context,
//                 const vertex_type& vertex,
//                 edge_type& edge) const {
//   not needed
//     }

void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    const vertex_type other = get_other_vertex(edge, vertex);

//		std::cout<<"scattering to "<<other.id()<<"\n";
if(other.data().type==14 && vertex.data().type == 14) {//category to category 
	our_msg for_cat;
	
	for(int i=0; i< temp1.size(); i++){
    if(temp1[i].from_last_art <12){
	   	triple tp2(temp1[i].vid,temp1[i].from_src +1,temp1[i].from_last_art+1);
      for_cat.collection.push_back(tp2); //updating distance
	   }
	 }
	if(for_cat.collection.size() !=0)
	context.signal(other,for_cat);
}
else if ( other.data().type == 0 && vertex.data().type==14){//category to article
    our_msg for_art;
		for(int i=0; i<temp1.size(); i++){
      if(other.data().vid_set.find( temp1[i].vid) != other.data().vid_set.end()){// It exists in the neighbours
       triple tp2(temp1[i].vid,temp1[i].from_src+1,temp1[i].from_last_art+1);
			 for_art.collection.push_back(tp2);
			}
		}
		if(for_art.collection.size() !=0)
		context.signal(other,for_art);
}
else if( other.data().type==14 && vertex.data().type ==0){//article to category
 our_msg for_cat1;
 triple tp2(vertex.id(),vertex.data().dist,0);
 for_cat1.collection.push_back(tp2);
 context.signal(other,for_cat1);
}

  } // end of scatter

// serialize
	 void save(graphlab::oarchive& oarc) const {
    oarc << temp1;
	   }
	
	// deserialize
	void load(graphlab::iarchive& iarc) {
	iarc >> temp1;
	}



};
 

/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
//------Output of the final graph----------
struct shortest_path_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
    std::stringstream strm;
    if(vtx.data().dist != std::numeric_limits<distance_type>::max() )
		strm << vtx.id() << "\t"<<"ns:"<<"\t"<<vtx.data().type<<"\t"<<vtx.data().dist<<"\n";
    /* boost::unordered_set<graphlab::vertex_id_type>::iterator it;
		 for( it = vtx.data().vid_set.begin(); it !=vtx.data().vid_set.end(); it++){
        strm<<*it<<" ";
		 }*/
		
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of shortest_path_writer

struct graph_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
   return "";
  }
  std::string save_edge(graph_type::edge_type e) { 
		 
		 std::stringstream strm;
		 strm<<e.source().data().type<<" "<<e.target().data().type<<"\n";
		 return strm.str();
		 }
}; // end of shortest_path_writer



struct max_deg_vertex_reducer: public graphlab::IS_POD_TYPE {
  size_t degree;
  graphlab::vertex_id_type vid;
  max_deg_vertex_reducer& operator+=(const max_deg_vertex_reducer& other) {
    if (degree < other.degree) {
      (*this) = other;
    }
    return (*this);
  }
};

max_deg_vertex_reducer find_max_deg_vertex(const graph_type::vertex_type vtx) {
  max_deg_vertex_reducer red;
  red.degree = vtx.num_in_edges() + vtx.num_out_edges();
  red.vid = vtx.id();
  return red;
}

bool line_parser_art(graph_type& graph, 
                 const std::string& filename, 
                 const std::string& textline) {
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
	graphlab::vertex_id_type tid;
  wiki_page_type type =1;  
 
	// first entry in the line is a vertex ID
  strm >> vid;
  strm >> tid;
	
 // std::cout<<textline<<std::endl;
	//dc.cout()<<textline<<std::endl;
  if(tid!=vid){

	//	std::cout<<"one \n"
 //	std::cout<<"adding edges \n";
		graph.add_edge(vid, tid);
    //std::cout<<"two \n";
	}
  return true;
}

bool line_parser_categ (graph_type& graph, 
                 const std::string& filename, 
                 const std::string& textline) {
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
	graphlab::vertex_id_type tid;
  wiki_page_type type =1;  
 
	// first entry in the line is a vertex ID
  strm >> vid;
  strm >> tid;
	
 // std::cout<<textline<<std::endl;
	//dc.cout()<<textline<<std::endl;
  if(tid!=vid){

	//	std::cout<<"one \n"
 //	std::cout<<"adding edges \n";
		graph.add_edge(vid, tid);
		graph.add_edge(tid,vid);
    //std::cout<<"two \n";
	}
  return true;
}
bool all_vertex_parser(graph_type& graph,
                  const std::string& filename,
                  const std::string& textline){
std::stringstream strm(textline);
graphlab::vertex_id_type vid;
wiki_page_type type;

// first entry in the line is a vertex ID
strm >> vid;
strm >> type;

// std::cout<<textline<<std::endl;
//dc.cout()<<textline<<std::endl;
  if(type==0 || type ==14)
  graph.add_vertex(vid, vertex_data(type));

	return true;

}
int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options 
    clopts("Single Source Shortest Path Algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  size_t powerlaw = 0;
  std::vector<graphlab::vertex_id_type> sources;
  bool max_degree_source = false;
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "graph format");
  clopts.attach_option("source", sources,
                       "The source vertices");
  clopts.attach_option("max_degree_source", max_degree_source,
                       "Add the vertex with maximum degree as a source");

  clopts.add_positional("source");

  clopts.attach_option("directed", DIRECTED_SSSP,
                       "Treat edges as directed.");

  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
 
  
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }


  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  
//	dc.cout() << "Loading graph in format: "<< format << std::endl;
  //Loading pages.txt for all the vertices!! 
  graph.load("hdfs://dsg2.crc.nd.edu/data/enwiki/page.txt" , all_vertex_parser);
  
  graph.load("hdfs://dsg2.crc.nd.edu/data/enwiki/pagelinks10m.txt", line_parser_art);
	
  graph.load("hdfs://dsg2.crc.nd.edu/data/enwiki/categorylinks.txt",line_parser_categ); 
	dc.cout()<<"Sources are---------after "<< sources[0]<<std::endl;
	// must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
            << "#edges:     " << graph.num_edges() << std::endl;



  if(sources.empty()) {
    if (max_degree_source == false) {
      dc.cout()
        << "No source vertex provided. Adding vertex 0 as source" 
        << std::endl;
      sources.push_back(0);
    }
  }

  if (max_degree_source) {
    max_deg_vertex_reducer v = graph.map_reduce_vertices<max_deg_vertex_reducer>(find_max_deg_vertex);
    dc.cout()
      << "No source vertex provided.  Using highest degree vertex " << v.vid << " as source."
      << std::endl;
    sources.push_back(v.vid);
  }
  // Running The Engine -------------------------------------------------------
  sources.clear();
  sources.push_back(303);	
	graphlab::omni_engine<add_neighbours> engine(dc, graph, exec_type, clopts);
  engine.signal_all();
	engine.start();
 /* 
	saveprefix = "/home/anambiar/algo_results/edges";
  if (saveprefix != "") {
    graph.save(saveprefix, graph_writer(),
               false,    // do not gzip
               false,     // save vertices
               true);   // do not save edges
  }
*/
	
	
	graphlab::omni_engine<main_algo> engine2(dc, graph, exec_type, clopts);
  our_msg init_msg;
 triple tp2(sources[0],0,0);
 init_msg.collection.push_back(tp2);
 
 
  engine2.signal(sources[0], init_msg);

  
  //Signal all the vertices in the source set
 /* for(size_t i = 0; i < sources.size(); ++i) {
    engine.signal(sources[i], min_distance_type(0));
  }
	*/
	//NO SIGNALLING CORRECT IT
  /*engine.signal(sources[0], min_distance_type(0));
	 *
	 * */
  engine2.start();
  const float runtime = engine2.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime
            << " seconds." << std::endl;



  // Save the final graph -----------------------------------------------------
  saveprefix = "/home/anambiar/algo_results/Yay";
  if (saveprefix != "") {
    graph.save(saveprefix, shortest_path_writer(),
               false,    // do not gzip
               true,     // save vertices
               false);   // do not save edges
  }
  
  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
	
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation



