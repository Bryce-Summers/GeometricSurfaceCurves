#include "halfEdgeMesh.h"

namespace CMU462 {

 bool Halfedge::isBoundary( void )
  // returns true if and only if this halfedge is on the boundary
  {
	return face()->isBoundary();
  }

  bool Edge::isBoundary( void )
  {
	return halfedge()->face()->isBoundary();
  }

  Vector3D Face::normal( void ) const
  {
     Vector3D N( 0., 0., 0. );

     HalfedgeCIter h = halfedge();
     do
     {
        Vector3D pi = h->vertex()->position;
        Vector3D pj = h->next()->vertex()->position;

        N += cross( pi, pj );

        h = h->next();
     }
     while( h != halfedge() );

     return N.unit();
  }

   void HalfedgeMesh :: build( const vector< vector<Index> >& polygons,
                               const vector<Vector3D>& vertexPositions )
   // This method initializes the halfedge data structure from a raw list of polygons,
   // where each input polygon is specified as a list of vertex indices.  The input
   // must describe a manifold, oriented surface, where the orientation of a polygon
   // is determined by the order of vertices in the list.  Polygons must have at least
   // three vertices.  Note that there are no special conditions on the vertex indices,
   // i.e., they do not have to start at 0 or 1, nor does the collection of indices have
   // to be contiguous.  Overall, this initializer is designed to be robust but perhaps
   // not incredibly fast (though of course this does not affect the performance of the
   // resulting data structure).  One could also implement faster initializers that
   // handle important special cases (e.g., all triangles, or data that is known to be
   // manifold).
   //
   // Since there are no strong conditions on the indices of polygons, we assume that
   // the list of vertex positions is given in lexicographic order (i.e., that the
   // lowest index appearing in any polygon corresponds to the first entry of the list
   // of positions and so on).
   {
      // define some types, to improve readability
      typedef vector<Index> IndexList;
      typedef IndexList::const_iterator IndexListCIter;
      typedef vector<IndexList> PolygonList;
      typedef PolygonList::const_iterator PolygonListCIter;
      typedef pair<Index,Index> IndexPair; // ordered pair of vertex indices, corresponding to an edge of an oriented polygon

      // Clear any existing elements.
       halfedges.clear();
        vertices.clear();
           edges.clear();
           faces.clear();
      boundaries.clear();

      // Since the vertices in our halfedge mesh are stored in a linked list,
      // we will temporarily need to keep track of the correspondence between
      // indices of vertices in our input and pointers to vertices in the new
      // mesh (which otherwise can't be accessed by index).  Note that since
      // we're using a general-purpose map (rather than, say, a vector), we can
      // be a bit more flexible about the indexing scheme: input vertex indices
      // aren't required to be 0-based or 1-based; in fact, the set of indices
      // doesn't even have to be contiguous.  Taking advantage of this fact makes
      // our conversion a bit more robust to different types of input, including
      // data that comes from a subset of a full mesh.
      map<Index,VertexIter> indexToVertex; // maps a vertex index to the corresponding vertex

      // Also store the vertex degree, i.e., the number of polygons that use each
      // vertex; this information will be used to check that the mesh is manifold.
      map<VertexIter,Size> vertexDegree;

      // First, we do some basic sanity checks on the input.
      for( PolygonListCIter p = polygons.begin(); p != polygons.end(); p++ )
      {
         if( p->size() < 3 )
         {
            // Refuse to build the mesh if any of the polygons have fewer than three vertices.
            // (Note that if we omit this check the code will still construct something fairly
            // meaningful for 1- and 2-point polygons, but enforcing this stricter requirement
            // on the input will help simplify code further downstream, since it can be certain
            // it doesn't have to check for these rather degenerate cases.)
            cerr << "Error converting polygons to halfedge mesh: each polygon must have at least three vertices." << endl;
            exit( 1 );
         }

         // We want to count the number of distinct vertex indices in this
         // polygon, to make sure it's the same as the number of vertices
         // in the polygon---if they disagree, then the polygon is not valid
         // (or at least, for simplicity we don't handle polygons of this type!).
         set<Index> polygonIndices;

         // loop over polygon vertices
         for( IndexListCIter i = p->begin(); i != p->end(); i++ )
         {
            polygonIndices.insert( *i );

            // allocate one vertex for each new index we encounter
            if( indexToVertex.find( *i ) == indexToVertex.end() )
            {
               VertexIter v = newVertex();
               v->halfedge() = halfedges.end(); // this vertex doesn't yet point to any halfedge
               indexToVertex[ *i ] = v;
               vertexDegree[ v ] = 1; // we've now seen this vertex only once
            }
            else
            {
               // keep track of the number of times we've seen this vertex
               vertexDegree[ indexToVertex[ *i ] ]++;
            }

         } // end loop over polygon vertices

         // check that all vertices of the current polygon are distinct
         Size degree = p->size(); // number of vertices in this polygon
         if( polygonIndices.size() < degree )
         {
            cerr << "Error converting polygons to halfedge mesh: one of the input polygons does not have distinct vertices!" << endl;
            cerr << "(vertex indices:";
            for( IndexListCIter i = p->begin(); i != p->end(); i++ )
            {
               cerr << " " << *i;
            }
            cerr << ")" << endl;
            exit( 1 );
         } // end check that polygon vertices are distinct

      } // end basic sanity checks on input

      // The number of vertices in the mesh is the
      // number of unique indices seen in the input.
      Size nVertices = indexToVertex.size();

      // The number of faces is just the number of polygons in the input.
      Size nFaces = polygons.size();
      faces.resize( nFaces ); // allocate storage for faces in our new mesh

      // We will store a map from ordered pairs of vertex indices to
      // the corresponding halfedge object in our new (halfedge) mesh;
      // this map gets constructed during the next loop over polygons.
      map<IndexPair,HalfedgeIter> pairToHalfedge;

      // Next, we actually build the halfedge connectivity by again looping over polygons
      PolygonListCIter p;
      FaceIter f;
      for( p = polygons.begin(), f = faces.begin();
           p != polygons.end();
           p++, f++ )
      {
         vector<HalfedgeIter> faceHalfedges; // cyclically ordered list of the half edges of this face
         Size degree = p->size(); // number of vertices in this polygon

         // loop over the halfedges of this face (equivalently, the ordered pairs of consecutive vertices)
         for( Index i = 0; i < degree; i++ )
         {
            Index a = (*p)[i]; // current index
            Index b = (*p)[(i+1)%degree]; // next index, in cyclic order
            IndexPair ab( a, b );
            HalfedgeIter hab;

            // check if this halfedge already exists; if so, we have a problem!
            if( pairToHalfedge.find( ab ) != pairToHalfedge.end() )
            {
               cerr << "Error converting polygons to halfedge mesh: found multiple oriented edges with indices (" << a << ", " << b << ")." << endl;
               cerr << "This means that either (i) more than two faces contain this edge (hence the surface is nonmanifold), or" << endl;
               cerr << "(ii) there are exactly two faces containing this edge, but they have the same orientation (hence the surface is" << endl;
               cerr << "not consistently oriented." << endl;
               exit( 1 );
            }
            else // otherwise, the halfedge hasn't been allocated yet
            {
               // so, we point this vertex pair to a new halfedge
               hab = newHalfedge();
               pairToHalfedge[ab] = hab;

               // link the new halfedge to its face
               hab->face() = f;
               hab->face()->halfedge() = hab;

               // also link it to its starting vertex
               hab->vertex() = indexToVertex[a];
               hab->vertex()->halfedge() = hab;

               // keep a list of halfedges in this face, so that we can later
               // link them together in a loop (via their "next" pointers)
               faceHalfedges.push_back( hab );
            }

            // Also, check if the twin of this halfedge has already been constructed (during
            // construction of a different face).  If so, link the twins together and allocate
            // their shared halfedge.  By the end of this pass over polygons, the only halfedges
            // that will not have a twin will hence be those that sit along the domain boundary.
            IndexPair ba( b, a );
            map<IndexPair,HalfedgeIter>::iterator iba = pairToHalfedge.find( ba );
            if( iba != pairToHalfedge.end() )
            {
               HalfedgeIter hba = iba->second;

               // link the twins
               hab->twin() = hba;
               hba->twin() = hab;

               // allocate and link their edge
               EdgeIter e = newEdge();
               hab->edge() = e;
               hba->edge() = e;
               e->halfedge() = hab;
            }
            else // If we didn't find a twin...
            {
               // ...mark this halfedge as being twinless by pointing
               // it to the end of the list of halfedges. If it remains
               // twinless by the end of the current loop over polygons,
               // it will be linked to a boundary face in the next pass.
               hab->twin() = halfedges.end();
            }

         } // end loop over the current polygon's halfedges

         // Now that all the halfedges of this face have been allocated,
         // we can link them together via their "next" pointers.
         for( Index i = 0; i < degree; i++ )
         {
            Index j = (i+1) % degree; // index of the next halfedge, in cyclic order
            faceHalfedges[i]->next() = faceHalfedges[j];
         }

      } // done building basic halfedge connectivity

      // For each vertex on the boundary, advance its halfedge pointer to one that is also on the boundary.
      for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
      {
         // loop over halfedges around vertex
         HalfedgeIter h = v->halfedge();
         do
         {
            if( h->twin() == halfedges.end() )
            {
               v->halfedge() = h;
               break;
            }

            h = h->twin()->next();
         }
         while( h != v->halfedge() ); // end loop over halfedges around vertex

      } // done advancing halfedge pointers for boundary vertices

      // Next we construct new faces for each boundary component.
      for( HalfedgeIter h = halfedgesBegin(); h != halfedgesEnd(); h++ ) // loop over all halfedges
      {
         // Any halfedge that does not yet have a twin is on the boundary of the domain.
         // If we follow the boundary around long enough we will of course eventually make a
         // closed loop; we can represent this boundary loop by a new face. To make clear the
         // distinction between faces and boundary loops, the boundary face will (i) have a flag
         // indicating that it is a boundary loop, and (ii) be stored in a list of boundaries,
         // rather than the usual list of faces.  The reason we need the both the flag *and* the
         // separate list is that faces are often accessed in two fundamentally different ways:
         // either by (i) local traversal of the neighborhood of some mesh element using the
         // halfedge structure, or (ii) global traversal of all faces (or boundary loops).
         if( h->twin() == halfedges.end() )
         {
            FaceIter b = newBoundary();
            vector<HalfedgeIter> boundaryHalfedges; // keep a list of halfedges along the boundary, so we can link them together

            // We now need to walk around the boundary, creating new
            // halfedges and edges along the boundary loop as we go.
            HalfedgeIter i = h;
            do
            {
               // create a twin, which becomes a halfedge of the boundary loop
               HalfedgeIter t = newHalfedge();
               boundaryHalfedges.push_back( t ); // keep a list of all boundary halfedges, in cyclic order
               i->twin() = t;
               t->twin() = i;
               t->face() = b;
               t->vertex() = i->next()->vertex();

               // create the shared edge
               EdgeIter e = newEdge();
               e->halfedge() = i;
               i->edge() = e;
               t->edge() = e;

               // Advance i to the next halfedge along the current boundary loop
               // by walking around its target vertex and stopping as soon as we
               // find a halfedge that does not yet have a twin defined.
               i = i->next();
               while( i != h && // we're done if we end up back at the beginning of the loop
                      i->twin() != halfedges.end() ) // otherwise, we're looking for the next twinless halfedge along the loop
               {
                  i = i->twin();
                  i = i->next();
               }
            }
            while( i != h );

            // The only pointers that still need to be set are the "next" pointers of the twins;
            // these we can set from the list of boundary halfedges, but we must use the opposite
            // order from the order in the list, since the orientation of the boundary loop is
            // opposite the orientation of the halfedges "inside" the domain boundary.
            Size degree = boundaryHalfedges.size();
            for( Index p = 0; p < degree; p++ )
            {
               Index q = (p-1+degree) % degree;
               boundaryHalfedges[p]->next() = boundaryHalfedges[q];
            }

         } // end construction of one of the boundary loops

         // Note that even though we are looping over all halfedges, we will still construct
         // the appropriate number of boundary loops (and not, say, one loop per boundary
         // halfedge).  The reason is that as we continue to iterate through halfedges, we
         // check whether their twin has been assigned, and since new twins may have been
         // assigned earlier in this loop, we will end up skipping many subsequent halfedges.

      } // done adding "virtual" faces corresponding to boundary loops

      // To make later traversal of the mesh easier, we will now advance the halfedge
      // associated with each vertex such that it refers to the *first* non-boundary
      // halfedge, rather than the last one.
      for( VertexIter v = verticesBegin(); v != verticesEnd(); v++ )
      {
         v->halfedge() = v->halfedge()->twin()->next();
      }

      // Finally, we check that all vertices are manifold.
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         // First check that this vertex is not a "floating" vertex;
         // if it is then we do not have a valid 2-manifold surface.
         if( v->halfedge() == halfedges.end() )
         {
            cerr << "Error converting polygons to halfedge mesh: some vertices are not referenced by any polygon." << endl;
            exit( 1 );
         }

         // Next, check that the number of halfedges emanating from this vertex in our half
         // edge data structure equals the number of polygons containing this vertex, which
         // we counted during our first pass over the mesh.  If not, then our vertex is not
         // a "fan" of polygons, but instead has some other (nonmanifold) structure.
         Size count = 0;
         HalfedgeIter h = v->halfedge();
         do
         {
            if( !h->face()->isBoundary() )
            {
               count++;
            }
            h = h->twin()->next();
         }
         while( h != v->halfedge() );

         if( count != vertexDegree[v] )
         {
            cerr << "Error converting polygons to halfedge mesh: at least one of the vertices is nonmanifold." << endl;
            exit( 1 );
         }
      } // end loop over vertices

      // Now that we have the connectivity, we copy the list of vertex
      // positions into member variables of the individual vertices.
      if( vertexPositions.size() != vertices.size() )
      {
         cerr << "Error converting polygons to halfedge mesh: number of vertex positions is different from the number of distinct vertices!" << endl;
         cerr << "(number of positions in input: " << vertexPositions.size() << ")" << endl;
         cerr << "(  number of vertices in mesh: " << vertices.size() << ")" << endl;
         exit( 1 );
      }
      // Since an STL map internally sorts its keys, we can iterate over the map from vertex indices to
      // vertex iterators to visit our (input) vertices in lexicographic order
      int i = 0;
      for( map<Index,VertexIter>::const_iterator e = indexToVertex.begin(); e != indexToVertex.end(); e++ )
      {
         // grab a pointer to the vertex associated with the current key (i.e., the current index)
         VertexIter v = e->second;

         // set the position of this vertex to the corresponding position in the input
         v->position = vertexPositions[ i ];
         i++;
      }

   } // end HalfedgeMesh::build()

   const HalfedgeMesh& HalfedgeMesh :: operator=( const HalfedgeMesh& mesh )
   // The assignment operator does a "deep" copy of the halfedge mesh data structure; in
   // other words, it makes new instances of each mesh element, and ensures that pointers
   // in the copy point to the newly allocated elements rather than elements in the original
   // mesh.  This behavior is especially important for making assignments, since the mesh
   // on the right-hand side of an assignment may be temporary (hence any pointers to elements
   // in this mesh will become invalid as soon as it is released.)
   {
      // Clear any existing elements.
       halfedges.clear();
        vertices.clear();
           edges.clear();
           faces.clear();
      boundaries.clear();

      // These maps will be used to identify elements of the old mesh
      // with elements of the new mesh.  (Note that we can use a single
      // map for both interior and boundary faces, because the map
      // doesn't care which list of faces these iterators come from.)
      map< HalfedgeCIter, HalfedgeIter > halfedgeOldToNew;
      map<   VertexCIter,   VertexIter >   vertexOldToNew;
      map<     EdgeCIter,     EdgeIter >     edgeOldToNew;
      map<     FaceCIter,     FaceIter >     faceOldToNew;

      // Copy geometry from the original mesh and create a map from
      // pointers in the original mesh to those in the new mesh.
      for( HalfedgeCIter h =      mesh.halfedgesBegin(); h !=  mesh.halfedgesEnd(); h++ ) halfedgeOldToNew[ h ] =  halfedges.insert(  halfedges.end(), *h );
      for(   VertexCIter v =       mesh.verticesBegin(); v !=   mesh.verticesEnd(); v++ )   vertexOldToNew[ v ] =   vertices.insert(   vertices.end(), *v );
      for(     EdgeCIter e =          mesh.edgesBegin(); e !=      mesh.edgesEnd(); e++ )     edgeOldToNew[ e ] =      edges.insert(      edges.end(), *e );
      for(     FaceCIter f =          mesh.facesBegin(); f !=      mesh.facesEnd(); f++ )     faceOldToNew[ f ] =      faces.insert(      faces.end(), *f );
      for(     FaceCIter b =     mesh.boundariesBegin(); b != mesh.boundariesEnd(); b++ )     faceOldToNew[ b ] = boundaries.insert( boundaries.end(), *b );

      // "Search and replace" old pointers with new ones.
      for( HalfedgeIter he = halfedgesBegin(); he != halfedgesEnd(); he++ )
      {
         he->next()   = halfedgeOldToNew[ he->next()   ];
         he->twin()   = halfedgeOldToNew[ he->twin()   ];
         he->vertex() =   vertexOldToNew[ he->vertex() ];
         he->edge()   =     edgeOldToNew[ he->edge()   ];
         he->face()   =     faceOldToNew[ he->face()   ];
      }
      for( VertexIter v =   verticesBegin(); v !=   verticesEnd(); v++ ) v->halfedge() = halfedgeOldToNew[ v->halfedge() ];
      for(   EdgeIter e =      edgesBegin(); e !=      edgesEnd(); e++ ) e->halfedge() = halfedgeOldToNew[ e->halfedge() ];
      for(   FaceIter f =      facesBegin(); f !=      facesEnd(); f++ ) f->halfedge() = halfedgeOldToNew[ f->halfedge() ];
      for(   FaceIter b = boundariesBegin(); b != boundariesEnd(); b++ ) b->halfedge() = halfedgeOldToNew[ b->halfedge() ];

      // Return a reference to the new mesh.
      return *this;
   }

   HalfedgeMesh :: HalfedgeMesh( const HalfedgeMesh& mesh )
   {
      *this = mesh;
   }


  // -- Local Operation implementations.


  
  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
   {

	 return doubleEdge(e0);

      // first handle boundary case
      if( e0->isBoundary() )
      {
         // get current elements
         HalfedgeIter h5 = e0->halfedge();
	 if( h5 -> isBoundary() ) h5 = h5 -> twin();
	 
         HalfedgeIter h4 = h5->next();
         HalfedgeIter h1 = h4->next();
         HalfedgeIter hb = h5->twin();
         VertexIter v1 = h1->vertex();
         VertexIter v3 = h5->vertex();
         VertexIter v4 = h4->vertex();
         EdgeIter e1 = h1->edge();
         EdgeIter e4 = h4->edge();
         FaceIter f1 = h1->face();
         FaceIter fb = hb->face();

         // get previous and next half edge along the boundary loop
         HalfedgeIter hbp = hb;
         do
         {
            hbp = hbp->twin()->next();
         }
         while( hbp->twin()->next() != hb );
         hbp = hbp->twin();
         HalfedgeIter hbn = hb->next();

         // allocate new elements
         VertexIter v5 = newVertex();
         EdgeIter e5 = newEdge();
         EdgeIter e7 = newEdge();
         EdgeIter e8 = newEdge();
         FaceIter f3 = newFace();
         FaceIter f6 = newFace();
         HalfedgeIter h8  = newHalfedge();
         HalfedgeIter h9  = newHalfedge();
         HalfedgeIter h11 = newHalfedge();
         HalfedgeIter h13 = newHalfedge();
         HalfedgeIter hb1 = newHalfedge();
         HalfedgeIter hb2 = newHalfedge();

         // set new vertex location
         v5->position = ( v3->position + v4->position ) / 2.;

         // connect new elements
         h1->setNeighbors( h11, h1->twin(), v1, e1, f3 );
         h4->setNeighbors( h9,  h4->twin(), v4, e4, f6 );
         h8->setNeighbors( h4,  hb1,      v5, e8, f6 );
         h9->setNeighbors( h8,  h13,      v1, e5, f6 );
         h11->setNeighbors( h13, hb2,     v3, e7, f3 );
         h13->setNeighbors( h1,  h9,      v5, e5, f3 );
         hb1->setNeighbors( hb2, h8,  v4, e8, fb );
         hb2->setNeighbors( hbn, h11, v5, e7, fb );
         hbp->next() = hb1;
         v1->halfedge() = h1;
         v3->halfedge() = h11;
         v4->halfedge() = h4;
         v5->halfedge() = h13;
         e1->halfedge() = h1;
         e4->halfedge() = h4;
         e5->halfedge() = h9;
         e7->halfedge() = h11;
         e8->halfedge() = h8;
         f3->halfedge() = h1;
         f6->halfedge() = h4;
         fb->halfedge() = hb1;

         // deallocate old elements
         deleteEdge( e0 );
         deleteFace( f1 );
         deleteHalfedge( h5 );
         deleteHalfedge( hb );

		 return v5;
      }
      else
      {
         // get current elements
         HalfedgeIter h5 = e0->halfedge();
         HalfedgeIter h6 = h5->twin();
         HalfedgeIter h4 = h5->next();
         HalfedgeIter h1 = h4->next();
         HalfedgeIter h3 = h6->next();
         HalfedgeIter h2 = h3->next();
         VertexIter v1 = h1->vertex();
         VertexIter v2 = h2->vertex();
         VertexIter v3 = h3->vertex();
         VertexIter v4 = h4->vertex();
         EdgeIter e1 = h1->edge();
         EdgeIter e2 = h2->edge();
         EdgeIter e3 = h3->edge();
         EdgeIter e4 = h4->edge();
         FaceIter f1 = h1->face();
         FaceIter f2 = h2->face();

         // allocate new elements
         VertexIter v5 = newVertex();
         EdgeIter e5 = newEdge();
         EdgeIter e6 = newEdge();
         EdgeIter e7 = newEdge();
         EdgeIter e8 = newEdge();
         FaceIter f3 = newFace();
         FaceIter f4 = newFace();
         FaceIter f5 = newFace();
         FaceIter f6 = newFace();
         HalfedgeIter h7  = newHalfedge();
         HalfedgeIter h8  = newHalfedge();
         HalfedgeIter h9  = newHalfedge();
         HalfedgeIter h10 = newHalfedge();
         HalfedgeIter h11 = newHalfedge();
         HalfedgeIter h12 = newHalfedge();
         HalfedgeIter h13 = newHalfedge();
         HalfedgeIter h14 = newHalfedge();

         // set new vertex location
         v5->position = ( v3->position + v4->position ) / 2.;

         // connect new elements
         h1->setNeighbors( h11, h1->twin(), v1, e1, f3 );
         h2->setNeighbors( h12, h2->twin(), v2, e2, f4 );
         h3->setNeighbors( h10, h3->twin(), v3, e3, f5 );
         h4->setNeighbors( h9,  h4->twin(), v4, e4, f6 );
         h7->setNeighbors( h3,  h11,      v5, e7, f5 );
         h8->setNeighbors( h4,  h12,      v5, e8, f6 );
         h9->setNeighbors( h8,  h13,      v1, e5, f6 );
         h10->setNeighbors( h7,  h14,      v2, e6, f5 );
         h11->setNeighbors( h13, h7,       v3, e7, f3 );
         h12->setNeighbors( h14, h8,       v4, e8, f4 );
         h13->setNeighbors( h1,  h9,       v5, e5, f3 );
         h14->setNeighbors( h2,  h10,      v5, e6, f4 );
         v1->halfedge() = h1;
         v2->halfedge() = h2;
         v3->halfedge() = h3;
         v4->halfedge() = h4;
         v5->halfedge() = h7;
         e1->halfedge() = h1;
         e2->halfedge() = h2;
         e3->halfedge() = h3;
         e4->halfedge() = h4;
         e5->halfedge() = h9;
         e6->halfedge() = h10;
         e7->halfedge() = h11;
         e8->halfedge() = h12;
         f3->halfedge() = h1;
         f4->halfedge() = h2;
         f5->halfedge() = h3;
         f6->halfedge() = h4;

         // deallocate old elements
         deleteEdge( e0 );
         deleteFace( f1 );
         deleteFace( f2 );
         deleteHalfedge( h5 );
         deleteHalfedge( h6 );
	 return v5;
      }
   }

   VertexIter HalfedgeMesh::collapseEdge( EdgeIter e )
   {

      HalfedgeIter he; // dummy iterator

      // handle boundary case first
      if( e->isBoundary() )
      {
         // get pointers to the original geometry
         HalfedgeIter h0 = e->halfedge(); if( h0->isBoundary() ) h0 = h0->twin();
         HalfedgeIter h1 = h0->next();
         HalfedgeIter h2 = h1->next();
         HalfedgeIter h1f = h1->twin();
         HalfedgeIter h2f = h2->twin();
         HalfedgeIter hb = h0->twin();
         HalfedgeIter hbn = hb->next();
         HalfedgeIter hbp = hb; do { hbp = hbp->twin()->next(); } while( hbp->twin()->next() != hb ); hbp = hbp->twin();
         VertexIter v0 = h2->vertex();
         VertexIter v1 = h0->vertex();
         VertexIter v2 = h1->vertex();
         EdgeIter e1 = h1->edge();
         EdgeIter e2 = h2->edge();
         FaceIter f = h0->face();
         FaceIter fb = hb->face();

         // create new mesh elements
         VertexIter v = newVertex();
         EdgeIter e0 = newEdge();

         // link together elements
         // vertices
         v0->halfedge() = h1f;
         v->halfedge() = hbn;
         v->position = ( v1->position + v2->position ) / 2.;

         // edges
         e0->halfedge() = h1f;
         // faces
         fb->halfedge() = hbn;
         // halfedges
         he = v1->halfedge(); do { he->vertex() = v; he = he->twin()->next(); } while( he != v1->halfedge() );
         he = v2->halfedge(); do { he->vertex() = v; he = he->twin()->next(); } while( he != v2->halfedge() );
         h1f->twin() = h2f;
         h2f->twin() = h1f;
         h1f->edge() = h2f->edge() = e0;
         hbp->next() = hbn;

         // remove old elements
         deleteVertex( v1 );
         deleteVertex( v2 );
         deleteEdge( e );
         deleteEdge( e1 );
         deleteEdge( e2 );
         deleteFace( f );
         deleteHalfedge( h0 );
         deleteHalfedge( h1 );
         deleteHalfedge( h2 );
         deleteHalfedge( hb );
		 return v;
      }
      else
      {
         // get fixed pointers to the original geometry
         HalfedgeIter h0 = e->halfedge();
         HalfedgeIter h1 = h0->twin();
         HalfedgeIter h00 = h0->next(), h01 = h00->next();
         HalfedgeIter h10 = h1->next(), h11 = h10->next();
         HalfedgeIter h00f = h00->twin();
         HalfedgeIter h01f = h01->twin();
         HalfedgeIter h10f = h10->twin();
         HalfedgeIter h11f = h11->twin();
         VertexIter v0 = h0->vertex();
         VertexIter v1 = h1->vertex();
         VertexIter v01 = h01->vertex();
         VertexIter v11 = h11->vertex();
         EdgeIter e00 = h00->edge();
         EdgeIter e01 = h01->edge();
         EdgeIter e10 = h10->edge();
         EdgeIter e11 = h11->edge();
         FaceIter f0 = h0->face();
         FaceIter f1 = h1->face();
         bool v0b = v0->isBoundary();
         bool v1b = v1->isBoundary();

         // check that the intersection of the 1-ring neighborhoods of the two
         // edge endpoints contains only the two vertices opposite the edge
         set<VertexCIter> neighbors;
         size_t n = 0;
         he = v0->halfedge(); do { n++; neighbors.insert( he->twin()->vertex() ); he = he->twin()->next(); } while( he != v0->halfedge() );
         he = v1->halfedge(); do { n++; neighbors.insert( he->twin()->vertex() ); he = he->twin()->next(); } while( he != v1->halfedge() );
         if( n - neighbors.size() != 2 )
         {
		   return verticesEnd();;
         }

         // create new mesh elements
         EdgeIter e0 = newEdge();
         EdgeIter e1 = newEdge();
         VertexIter v = newVertex();

         // link old elements to new ones
         // vertices
         v->halfedge() = h01f;
         // edges
         e0->halfedge() = h01f;
         e1->halfedge() = h11f;
         // faces
         // (nothing to do)
         // halfedges
         he = v0->halfedge(); do { he->vertex() = v; he = he->twin()->next(); } while( he != v0->halfedge() );
         he = v1->halfedge(); do { he->vertex() = v; he = he->twin()->next(); } while( he != v1->halfedge() );
         v01->halfedge() = h00f;
         v11->halfedge() = h10f;
         h00f->edge() = e0;
         h01f->edge() = e0;
         h10f->edge() = e1;
         h11f->edge() = e1;
         h00f->twin() = h01f;
         h01f->twin() = h00f;
         h10f->twin() = h11f;
         h11f->twin() = h10f;

         // if both vertices were on the bounday, put new vertex
         // at average of neighbors; otherwise, put it at the boundary
         if( v0b )
         {
            v->position = v0->position;
         }
         else if( v1b )
         {
            v->position = v1->position;
         }
         else
         {
            v->position = (v0->position+v1->position)/2.;
         }

         // // put new vertex at centroid of neighbors
         // double k = 0.;
         // v->position = Vector( 0., 0., 0. );
         // he = v->halfedge();
         // do
         // {
         //    k += 1.;
         //    v->position += he->twin()->vertex()->position;
         //    he = he->twin()->next();
         // }
         // while( he != v->halfedge() );
         // v->position /= k;

         // remove old mesh elements
         deleteVertex( v0 );
         deleteVertex( v1 );
         deleteEdge( e );
         deleteEdge( e00 );
         deleteEdge( e01 );
         deleteEdge( e10 );
         deleteEdge( e11 );
         deleteFace( f0 );
         deleteFace( f1 );
         deleteHalfedge( h0 );
         deleteHalfedge( h00 );
         deleteHalfedge( h01 );
         deleteHalfedge( h1 );
         deleteHalfedge( h10 );
         deleteHalfedge( h11 );

		 return v;
      }

   }

   EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
   {

      // skip boundary edges
      if( e0->isBoundary() ) return edgesEnd();

      // get current elements
      HalfedgeIter h1 = e0->halfedge();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h2->next();
      HalfedgeIter h4 = h1->twin();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h5->next();
      VertexIter v0 = h4->vertex();
      VertexIter v1 = h1->vertex();
      VertexIter v2 = h3->vertex();
      VertexIter v3 = h6->vertex();
      FaceIter f0 = h1->face();
      FaceIter f1 = h4->face();

      // recconnect current elements
      h1->setNeighbors( h6, h4, v2, e0, f0 );
      h2->setNeighbors( h1, h2->twin(), v0, h2->edge(), f0 );
      h3->setNeighbors( h5, h3->twin(), v2, h3->edge(), f1 );
      h4->setNeighbors( h3, h1, v3, e0, f1 );
      h5->setNeighbors( h4, h5->twin(), v1, h5->edge(), f1 );
      h6->setNeighbors( h2, h6->twin(), v3, h6->edge(), f0 );
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h1;
      v3->halfedge() = h4;
      e0->halfedge() = h1;
      f0->halfedge() = h1;
      f1->halfedge() = h4;

	  return e0;
   }

  // Transforms 1 edge into 2. Returns a pointer to the middle vertex.
  // Needs to make sure that the two new edges are labeled as new.
  VertexIter HalfedgeMesh::doubleEdge(EdgeIter e1)
  {

	// Find Relevant Half Edges.
	HalfedgeIter h1 = e1 -> halfedge();
	HalfedgeIter h2 = h1 -> twin();

	// Allocate New Elements.
	VertexIter v = newVertex();
	HalfedgeIter h1_next = newHalfedge();
	HalfedgeIter h2_next = newHalfedge();
	EdgeIter e2 = newEdge();

	e1 -> isNew = true;
	e2 -> isNew = true;

	// 6 Linkings on one side.
	h1_next -> next() = h1 -> next();
	h1 -> next() = h1_next;
	h1_next -> face() = h1 -> face();
	h1_next -> vertex() = v;
	h1_next -> twin() = h2;
	h1 -> twin() = h2_next;

	// 6 linkings on the other.
	h2_next -> next() = h2 -> next();
	h2 -> next() = h2_next;
	h2_next -> face() = h2 -> face();
	h2_next -> vertex() = v;
	h2 -> twin() = h1 -> next();
	h2_next -> twin() = h1;

	e1 -> halfedge() = h1;
	e2 -> halfedge() = h2;
	h1 -> edge() = e1;
	h2 -> edge() = e2;
	h1_next -> edge() = e2;
	h2_next -> edge() = e1;

	Vector3D pos1 = h1 -> vertex() -> position;
	Vector3D pos2 = h2 -> vertex() -> position;
	v -> position = (pos1 + pos2)/2;
	v -> halfedge() = h1_next;

	return v;
  }


  

} // End of CMU 462 namespace.
