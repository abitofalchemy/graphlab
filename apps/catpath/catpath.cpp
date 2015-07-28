/** S. Aguinaga
 ** ** ** ** ** ** */

/*-------*---------*---------*---------*---------*---------*---------*---------*
 * Author: s. aguinaga saguinag (at) nd dot edu
 * * This work is derive from:
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
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

bool DIRECTED_SSSP = true;

/**
 *  * \brief The type used to measure distances in the graph.
 *   */
typedef float distance_type;

struct Quad {
	graphlab::vertex_id_type dest_art;
	distance_type from_last_art;
	distance_type from_src;
	//  graphlab::vertex_id_type last_cat;
	//	graphlab::vertex_id_type last_art;
	Quad(): dest_art(0) {}

	Quad(graphlab::vertex_id_type a, distance_type b, distance_type c){
		dest_art =a;
		from_last_art = c;
		from_src = b;
	//	last_art = _last_art;
  	//  last_cat = d;
	}

	void save(graphlab::oarchive& oarc) const { oarc << dest_art << from_last_art << from_src;}
	void load(graphlab::iarchive& iarc) { iarc >> dest_art >> from_last_art >> from_src;}

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
struct vertex_data : graphlab::IS_POD_TYPE {
    distance_type dist;
	vertex_data(distance_type dist = std::numeric_limits<distance_type>::max()) :
            dist(dist) { }
}; // end of vertex data



/**
 *  * \brief The distance associated with the edge.
 *   */
struct edge_data : graphlab::IS_POD_TYPE {
	distance_type dist;
	edge_data(distance_type dist = 1) : dist(dist) { }
}; // end of edge data


/**
 *  * \brief The graph type encodes the distances between vertices and
 *   * edges
 *    */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

inline graph_type::vertex_type get_other_vertex(const graph_type::edge_type& edge,
                     const graph_type::vertex_type& vertex) {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
}

/** 
 ** main_algo (catpath? from zliu8)
 **
 **/
// class main_algo:
// 	public graphlab::ivertex_program<graph_type,
// 	graphlab::empty, 
// 	our_msg >  // our_msg is my message type not gather type 
// 
// {
// 
// 	std::vector<Quad> temp1;
// 	public:
// 	void init(icontext_type& context, const vertex_type& vertex,
// 			const our_msg& msg) {
// 		temp1 = msg.collection;
// 
// 	} 
// 
// 
// 	// Gather on all edges
// 	edge_dir_type gather_edges(icontext_type& context,
// 			const vertex_type& vertex) const {
// 		return graphlab::NO_EDGES;
// 	}
// 
// 	void apply(icontext_type& context, vertex_type& vertex,
// 			const graphlab::empty& empty) 
// 	{
// 
// 		distance_type tp = std::numeric_limits<distance_type>::max() ;
// 
// 		if(vertex.data().sent == true){
// 			vertex.data().isDead = true;
// 		} else if(vertex.data().type==0){ // if a article
// 
// 			//std::cout<<"Applying on "<<vertex.id()<< " --current "<<vertex.data().dist<<std::endl;
// 			for(int i=0; i<temp1.size();i++){
// 				//gets minimum distance from incoming category edges
// 				if(temp1[i].from_src < tp){
// 					tp=temp1[i].from_src;
// 				}
// 			}
// 			//is the minimum we just found, less than the distance we currently store? If so, store new minimum.
// 			if( tp < vertex.data().dist){
// 				vertex.data().dist = tp;
// 				//		std::cout<<"Setting distance of "<<vertex.id()<<"-"<<tp<<std::endl;
// 				vertex.data().sent = true;
// 			}
// 
// 
// 			//end of if an article	
// 		} else if(vertex.data().type==14){//its a category
// 			// dont do anything
// 			
// 			std::map<graphlab::vertex_id_type, Triple> temp;
// 								
// 			for(unsigned int i=0; i < temp1.size(); i++)
// 			{
// 				vertex.data().seen.insert(temp1[i].dest_art);
// 
// 				if(vertex.data().msg_q.size() > 0){
// 					vertex.data().msg_q.pop_front(); //removed in the previous scatter
// 				}
// 				
// 				vertex.data().msg_q.push_back( temp1[i] ) ;
// 
// 			}
// 			
// 		}
// 	} // end of apply
// 
// 	edge_dir_type scatter_edges(icontext_type& context,
// 			const vertex_type& vertex) const {
// 		return graphlab::OUT_EDGES;
// 	}
// 
// 	// void scatter(icontext_type& context,
// 	//                 const vertex_type& vertex,
// 	//                 edge_type& edge) const {
// 	//   not needed
// 	//     }
// 
// 	void scatter(icontext_type& context, const vertex_type& vertex,
// 			edge_type& edge) const 
// 	{
// 		const vertex_type other = get_other_vertex(edge, vertex);
// 		//		std::cout<<"scattering to "<<other.id()<<"\n";
// 		if(other.data().type==14 && vertex.data().type == 14) {//category to category 
// 			our_msg for_cat;
// 
// 			for(int i=0; i < 1/*temp1.size()*/; i++)
// 			{
// 				//if(temp1[i].from_last_art < 12 && temp1[i].last_cat != other.id() )
// 				{
// 
// 					Quad msg = vertex.data().msg_q.front() ;
// 					if(other.data().seen.find(msg.dest_art) != other.data().seen.end())
// 					{
// 						msg.from_src += 1;
// 						msg.from_last_art += 1;
// 
// 						for_cat.collection.push_back(msg); //updating distance
// 					}
// 				}
// 			}
// 		 	if(for_cat.collection.size() !=0){
// 				context.signal(other, for_cat);
// 			}
//    
// 		} else if ( other.data().type == 0 && vertex.data().type==14 && other.data().isDead == false){//category to article
//    
// 		  our_msg for_art;
// 			for(int i=0; i<temp1.size(); i++){
// 				//If it exists in the add_neighbours 
// 				if(other.data().vid_set.find( temp1[i].dest_art ) != other.data().vid_set.end()// && 
// 					//(other.data().seen.find(temp1[i].dest_art) != other.data().seen.end()) )
// 					)
// 				{			
// 					Quad tp2(temp1[i].dest_art, temp1[i].from_src+1, temp1[i].from_last_art+1 );
// 					if(temp1[i].from_src != std::numeric_limits<distance_type>::max())
// 					{
// 						for_art.collection.push_back(tp2);
// 					}
// 				}
// 			}
// 			if(for_art.collection.size() !=0){
// 				context.signal(other,for_art);
// 			}
// 		} else if( other.data().type==14 && vertex.data().type ==0 && vertex.data().isDead == false){//article to category
// 		  //std::cout<<"Required \n";
// 
// 			our_msg for_cat1;
// 		  Quad tp2(other.id(), vertex.data().dist, 0);
// 		  for_cat1.collection.push_back(tp2);
// 		  if(vertex.data().dist != std::numeric_limits<distance_type>::max()) 
// 			{
// 			  context.signal(other, for_cat1);
// 			}
// 		}
// 	
// 	} // end of scatter
// 
// // serialize
// 	 void save(graphlab::oarchive& oarc) const {
//     oarc << temp1;
// 	   }
// 	
// 	// deserialize
// 	void load(graphlab::iarchive& iarc) {
// 		iarc >> temp1;
// 	}
// 
// 
// 
// }; // ends main algo


// Class that used with the omni engine
//class spaths :
// Shortest Paths: Derived from bshi 
// HDTM/src/graphlab/apps/hdtm_analysis/hdtm_ana.cpp
class shortest_paths:
	public graphlab::ivertex_program</*cate_*/graph_type, 
									 graphlab::empty, 
									 our_msg/*sp_message_type*/>
{
	std::vector<Quad> temp1;
	public:
	void init(icontext_type& context, const vertex_type& vertex, const our_msg& msg) 
	{
		temp1 = msg.collection;
		}
	
	// Gather on all edges:
	edge_dir_type gather_edges(icontext_type& 		context,
							   const vertex_type&   vertex) const 
	{	return graphlab:NO_EDGES; }
	
	void apply(	icontext_type&  context, 
				vertex_type&    vertex,
				const graphlab::empty& empty) 
	{
		distance_type tp = std::numeric_limits<distance_type>::max(); // init to inf?

		if(vertex.data().sent == true){
			vertex.data().isDead = true;
		} else if(vertex.data().type==0){ // if a article, but is this might not be avail
			//std::cout<<"Applying on "<<vertex.id()<< " --current "<<vertex.data().dist<<std::endl;
			std::cout << "temp1.size() " << temp1.size() << std::endl;
			for(int i=0; i<temp1.size();i++){
				//gets minimum distance from incoming category edges
				if(temp1[i].from_src < tp){
					tp=temp1[i].from_src;
				}
			}
			//is the minimum we just found, less than the distance we currently store? If so, store new minimum.
			if( tp < vertex.data().dist){
				vertex.data().dist = tp;
				//		std::cout<<"Setting distance of "<<vertex.id()<<"-"<<tp<<std::endl;
				vertex.data().sent = true;
			}


			//end of if an article	
		} else if(vertex.data().type==14){//its a category
			// dont do anything
			
			std::map<graphlab::vertex_id_type, Triple> temp;
								
			for(unsigned int i=0; i < temp1.size(); i++)
			{
				vertex.data().seen.insert(temp1[i].dest_art);

				if(vertex.data().msg_q.size() > 0){
					vertex.data().msg_q.pop_front(); //removed in the previous scatter
				}
				
				vertex.data().msg_q.push_back( temp1[i] ) ;

			}
			
		}
	} // end of apply
		edge_dir_type scatter_edges(icontext_type& context,
			const vertex_type& vertex) const {
		return graphlab::OUT_EDGES;
	}

	// Scatter method
	void scatter(icontext_type& context, const vertex_type& vertex,
			edge_type& edge) const 
	{
		const vertex_type other = get_other_vertex(edge, vertex);
		//		std::cout<<"scattering to "<<other.id()<<"\n";
		if(other.data().type==14 && vertex.data().type == 14) {//category to category 
			our_msg for_cat;

			for(int i=0; i < 1/*temp1.size()*/; i++)
			{
				//if(temp1[i].from_last_art < 12 && temp1[i].last_cat != other.id() )
				{

					Quad msg = vertex.data().msg_q.front() ;
					if(other.data().seen.find(msg.dest_art) != other.data().seen.end())
					{
						msg.from_src += 1;
						msg.from_last_art += 1;

						for_cat.collection.push_back(msg); //updating distance
					}
				}
			}
		 	if(for_cat.collection.size() !=0){
				context.signal(other, for_cat);
			}
   
		} else if ( other.data().type == 0 && vertex.data().type==14 && other.data().isDead == false){//category to article
   
		  our_msg for_art;
			for(int i=0; i<temp1.size(); i++){
				//If it exists in the add_neighbours 
				if(other.data().vid_set.find( temp1[i].dest_art ) != other.data().vid_set.end()// && 
					//(other.data().seen.find(temp1[i].dest_art) != other.data().seen.end()) )
					)
				{			
					Quad tp2(temp1[i].dest_art, temp1[i].from_src+1, temp1[i].from_last_art+1 );
					if(temp1[i].from_src != std::numeric_limits<distance_type>::max())
					{
						for_art.collection.push_back(tp2);
					}
				}
			}
			if(for_art.collection.size() !=0){
				context.signal(other,for_art);
			}
		} else if( other.data().type==14 && vertex.data().type ==0 && vertex.data().isDead == false){//article to category
			//std::cout<<"Required \n";
			our_msg for_cat1;
			Quad tp2(other.id(), vertex.data().dist, 0);
			for_cat1.collection.push_back(tp2);
			if(vertex.data().dist != std::numeric_limits<distance_type>::max()) 
			{
				context.signal(other, for_cat1);
				}
		}
	
	} // end of scatter

	// serialize
	void save(graphlab::oarchive& oarc) const {	oarc << temp1;}
	
	// deserialize
	void load(graphlab::iarchive& iarc) { iarc >> temp1;}



}; // ends main algo
	


int main( int argc, char** argv) { 
  // Initialization block 
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options
  graphlab::command_line_options 
        clopts("Catpath: Augmented Shortest Path Algorithm.");
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
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant catpath to a "
                       "sequence of files with prefix saveprefix");
  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  // Begin
  // Build the graph
  graph_type graph(dc, clopts);
  if (powerlaw >0 ) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2, 100000000);
    } else if (graph_dir.length() > 0) { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load_format(graph_dir, format);
    } else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
    }
  
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
            << "#edges:     " << graph.num_edges() << std::endl;
  
  // Source nodes
  if(sources.empty()) {
    if (max_degree_source == false) {
      dc.cout()
      << "No source vertex provided. Adding vertex 0 as source" 
      << std::endl;
      sources.push_back(0);
      }
    }
  
  // Running The Engine
  graphlab::omni_engine<shortest_paths> engine(dc, graph, exec_type, clopts);



  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS; 
} // end of main
