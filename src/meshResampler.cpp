/*
 * Student solution for CMU 15-462 Project 2 (MeshEdit)
 *
 * Implemented by ____ on ____.
 *
 */

#include "meshResampler.h"
#include "mutablePriorityQueue.h"

namespace CMU462
{

   void MeshResampler::upsample( HalfedgeMesh& mesh )
   {
	 // -- Catmull - Clark Subdivision.
	 /* Step 1. Compute a 'face point' for every face. Store it in the face.
	  *         Set its position to be the weighted average of all positions on the face.
	  *
	  * Step 2. Split every edge in the mesh into 2.
	  *         Set its position to be be the average of the associated two face points and original edge endpoints.
	  *         While we are at it, we might as well add our midpoint contribution 'R' that will be used in Step 4.
	  *
	  * Step 3. Connect every face point to the set of neighboring edge points.
	  *         Two Passes.
	  *
	  * Step 4. Compute new positions for each of the original mesh points (P).
	  * F = average of neighboring face points. (n of them).
	  * R = average of all edge_points for edges touching P.
	  * new position = [ F + 2R + (n - 3)P ] / n
	  *
	  */


	 // Set all edges to 'not new'.
	 for(EdgeIter edge = mesh.edgesBegin(); edge != mesh.edgesEnd(); edge++)
	 {
	   edge -> isNew = false;
	 }

	 // Set all vertices to 'not new'.
	 for(VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++)
	 {
	   vertex -> isNew = false;
	   vertex -> newPosition = Vector3D(0, 0, 0);
	   vertex -> one_neighborhood_midpoint_sum = Vector3D(0, 0, 0);
	 }

	 // Set all faces to 'not new'.
	 for(FaceIter face = mesh.facesBegin(); face != mesh.facesEnd(); face++)
	 {
	   face -> isNew = false;
	 }

	 cout << "Step 1!";

	 // -- Step 1. Face points.
	 for(FaceIter face = mesh.facesBegin(); face != mesh.facesEnd(); face++)
	 {
	     Vector3D & face_point = face -> face_point;
		 face_point = Vector3D(0, 0, 0);

		 // Compute weighted average of vertices on the face.
		 int n = 0;
		 HalfedgeIter h0 = face -> halfedge();
		 HalfedgeIter h = h0;

		 do
		 {
		   face_point += h -> vertex() -> position;
		   n++;

		   // Go to the next halfedge.
		   h = h->next();
		 }
		 while(h != h0);

		 face_point /= n;
	 }

	 cout << "Step 2!" << endl;

	 // -- Step 2. Edge points.

	 EdgeIter edge_iter = mesh.edgesBegin();
	 while(edge_iter != mesh.edgesEnd())
	 {

	   // Don't split the new edges.
	   if(edge_iter -> isNew)
	   {
		 edge_iter++;
		 continue;
	   }

	   EdgeIter edge_iter_next = edge_iter;
	   edge_iter_next++;

	   HalfedgeIter h1 = edge_iter -> halfedge();
	   HalfedgeIter h2 = h1 -> twin();

	   VertexIter v1 = h1 -> vertex();
	   VertexIter v2 = h2 -> vertex();

	   Vector3D endpoint1 = v1 -> position;
	   Vector3D endpoint2 = v2 -> position;
	   Vector3D face_point1 = h1 -> face() -> face_point;
	   Vector3D face_point2 = h2 -> face() -> face_point;

	   // Mutate the halfedge mesh by doubling the edges through midpoint splitting.
	   VertexIter edge_point   = mesh.doubleEdge(edge_iter);

	   edge_point -> isNew = true;


	   // Update the edge point to be the average of the original edge midpoint and
	   // the midpoint of the line segment between the neighboring face points.
	   // ((E1 + E1)/2 + (F1 + F2)/2)/2 = (E1 + E2 + F1 + F2)/4
	   edge_point -> position  = endpoint1 + endpoint2 + face_point1 + face_point2;
	   edge_point -> position /= 4.0;


	   // Add the edge point contribution to each original vertex.
	   v1 -> one_neighborhood_midpoint_sum  +=  edge_point -> position;
	   v2 -> one_neighborhood_midpoint_sum  +=  edge_point -> position;


	   // Iterate to the next edge.
	   edge_iter = edge_iter_next;
	 }

	 cout << "Step 3!" << endl;


	 // -- Step 3. Connect the Face points to the edge points.
	 FaceIter face_iter = mesh.facesBegin();
	 while(face_iter != mesh.facesEnd())
	 {
	   // We only need to link up old faces.
	   if(face_iter -> isNew)
	   {
		 face_iter++;
		 continue;
	   }

	   FaceIter face_iter_next = face_iter;
	   face_iter_next++;

	   // Allocate the new vertex for the face.
       VertexIter face_vertex  = mesh.newVertex();
	   face_vertex -> position = face_iter -> face_point;
	   face_vertex -> isNew = true;

	   // -- Go around the face, linking up the new 'edge points' to the face point.
	   HalfedgeIter h0 = face_iter -> halfedge();

	   // We want to ensure that we are on the right parity.
	   // the first half edge should be pointing towards a new edge point,
	   // instead of away from it so that we are connecting the correct edges.
	   if(h0 -> vertex() -> isNew)
	   {
		 h0++;
	   }

	   HalfedgeIter edge_current = h0;


	   // -- 3A Pass 1. Link every edge_next pointer in the loop to a new face.
	   //               Also do some initial linking for the new edge and halfedges.

	   do
	   {
		 HalfedgeIter edge_next = edge_current -> next();

		 // -- Allocate new elements.
		 FaceIter new_face = mesh.newFace();
		 new_face -> isNew = true;
		 HalfedgeIter h1 =  mesh.newHalfedge();
		 HalfedgeIter h2 =  mesh.newHalfedge();
		 EdgeIter e1 = mesh.newEdge();

		 // Assignment links.
		 edge_current -> next() = h1;
		 h1 -> vertex() = edge_next -> vertex();
		 h1 -> twin() = h2;
		 e1 -> halfedge() = h1;
		 h1 -> edge() = e1;
		 h2 -> twin() = h1;
		 // Note: 'face_vertex' is the location of the face point.
		 h2 -> vertex() = face_vertex;
		 h2 -> face()   = new_face;
		 h2 -> next()   = edge_next;
		 edge_next -> face() = new_face;
		 new_face -> halfedge() = edge_next;


		 // Add the Face point contributions to every original vertex.
		 edge_current->vertex() -> newPosition += face_vertex -> position;


		 // Go on to the next set of 2 edges,
		 // since we are only connecting the new edge points.
		 edge_current = edge_next -> next();

	   }while(edge_current != h0);


	   // Give the center face vertex a link to one of its outgoing half edges.
	   face_vertex -> halfedge() = edge_current -> next() -> twin();


	   // Step 3B. Complete the linking by linking the looping circle of halfedges.
	   // ASSERTION(edge_current == h0);
	   do
	   {
		 // Find Relevant elements.
		 HalfedgeIter h1 = edge_current -> next();
		 HalfedgeIter h2 = h1 -> twin();
		 HalfedgeIter edge_next = h2 -> next();
		 FaceIter f1 = edge_next -> face();
		 HalfedgeIter next_edge_current = edge_next -> next();
		 HalfedgeIter next_h1 = next_edge_current   -> next();

		 next_edge_current -> face() = f1;
		 next_h1 -> face() = f1;
		 next_h1 -> next() = h2;

		 // Go on to the next set of 2 edges,
		 // since we are only connecting the new edge points.
		 edge_current = next_edge_current;

	   }while(edge_current != h0);

	   // Delete the original face.
	   mesh.deleteFace(face_iter);

	   // Iterate to the next face.
	   face_iter = face_iter_next;
	 }

	 cout << "Step 4!\n";

	 /* Step 4. Do a whole bunch of averaging.
	  * Compute new positions for each of the original mesh points (P).
	  * F = average of neighboring face points. (n of them)
	  * R = average of all edge_points for edges touching P (n of them).
	  *      these edge points were calculated during Step 2.
	  * new position = [ F + 2R + (n - 3)P ] / n
	  */
	 for(VertexIter vertex = mesh.verticesBegin(); vertex != mesh.verticesEnd(); vertex++)
	 {
	   // Don't update the positions of new points.
	   if(vertex -> isNew)
	   {
		 continue;
	   }

	   // The F value is computed during step 2.
	   Vector3D F = vertex -> newPosition;
	   Vector3D R = vertex -> one_neighborhood_midpoint_sum;
	   Vector3D P = vertex -> position;
	   int n = vertex -> degree();

	   // We have summed up the contributions at disparate points in time.
	   // We therefore need to divide by the size of the one neighborhood to convert these
	   // values to averages.
	   F /= n;
	   R /= n;

	   // Update the vertex position as an affine combination of F, R and P.
	   vertex -> position = (F + 2*R + (n - 3)*P)/n;

	 }
	 //ZZZ
	 cout << "Catmull - Clark Subdivision is complete!\n";

   }// -- End of Catmull - Clark Subdivision.

   // Given an edge, the constructor for EdgeRecord finds the
   // optimal point associated with the edge's current quadric,
   // and assigns this edge a cost based on how much quadric
   // error is observed at this optimal point.
   EdgeRecord::EdgeRecord( EdgeIter& _edge )
   : edge( _edge )
   {
      // TODO Compute the combined quadric from the edge endpoints.


      // TODO Build the 3x3 linear system whose solution minimizes
      // the quadric error associated with these two endpoints.


      // TODO Use this system to solve for the optimal position, and
      // TODO store it in EdgeRecord::optimalPoint.


      // TODO Also store the cost associated with collapsing this edge
      // TODO in EdgeRecord::Cost.

   }

   void MeshResampler::downsample( HalfedgeMesh& mesh )
   {
      // TODO Compute initial quadrics for each face by simply writing the plane
      // equation for the face in homogeneous coordinates.  These quadrics should
      // be stored in Face::quadric


      // TODO Compute an initial quadric for each vertex as the sum of the quadrics
      // associated with the incident faces, storing it in Vertex::quadric


      // TODO Build a priority queue of edges according to their quadric error cost,
      // TODO i.e., by building an EdgeRecord for each edge and sticking it in the queue.


      // TODO Until we reach the target edge budget, collapse the best edge.  Remember
      // TODO to remove from the queue any edge that touches the collapsing edge BEFORE
      // TODO it gets collapsed, and add back into the queue any edge touching the collapsed
      // TODO vertex AFTER it's been collapsed.  Also remember to assign a quadric to the
      // TODO collapsed vertex, and to pop the collapsed edge off the top of the queue.
   }

   void Vertex::computeCentroid( void )
   {
      // TODO Compute the average position of all neighbors of this vertex, and
      // TODO store it in Vertex::centroid.  This value will be used for resampling.
   }

   Vector3D Vertex::normal( void ) const
   // TODO Returns an approximate unit normal at this vertex, computed by
   // TODO taking the area-weighted average of the normals of neighboring
   // TODO triangles, then normalizing.
   {
      // TODO Compute and return the area-weighted unit normal.
      return Vector3D();
   }

   void MeshResampler::resample( HalfedgeMesh& mesh )
   {
      const int nIterations = 5;
      const int nSmoothingIterations = 20;
      EdgeIter e;

      // Compute the mean edge length; this will be the target length for remeshing.
      double meanEdgeLength = 0.;
      for( e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
      {
         meanEdgeLength += e->length();
      }
      meanEdgeLength /= (double) mesh.nEdges();
      meanEdgeLength *= .9;
      const double maxLength = meanEdgeLength * 4./3.;
      const double minLength = meanEdgeLength * 4./5.;

      for( int iteration = 0; iteration < nIterations; iteration++ )
      {
         cerr << "Adjusting edge lengths..." << endl;

         // We first try to get more uniform edge lengths by splitting edges longer
         // than the threshold and collapsing edges shorter than the threshold.
         e = mesh.edgesBegin();
         while( e != mesh.edgesEnd() )
         {
            double length = e->length();

            // After splitting or collapsing, the current edge may no longer
            // exist; therefore, we need to grab a pointer to the next edge NOW.
            EdgeIter nextEdge = e;
            nextEdge++;

            if( length > maxLength ) // Split edges that are longer than the target length.
            {
               mesh.splitEdge( e );
            }

            e = nextEdge;
         }

         // We first try to get more uniform edge lengths by splitting edges longer
         // than the threshold and collapsing edges shorter than the threshold.
         e = mesh.edgesBegin();
         while( e != mesh.edgesEnd() )
         {
            double length = e->length();

            // After splitting or collapsing, the current edge may no longer
            // exist; therefore, we need to grab a pointer to the next edge NOW.
            EdgeIter nextEdge = e;
            nextEdge++;

            if( length < minLength ) // Collapse edges that are shorter than the target length.
            {
               // Collapsing an edge may invalidate not only this edge, but
               // also any of the four edges in its immediate neighborhood.
               // Therefore, we will advance the next edge to the first edge
               // that is NOT one of these edges.
               EdgeIter e1 = e->halfedge()->next()->edge();
               EdgeIter e2 = e->halfedge()->next()->next()->edge();
               EdgeIter e3 = e->halfedge()->twin()->next()->edge();
               EdgeIter e4 = e->halfedge()->twin()->next()->next()->edge();
               while( nextEdge == e  ||
                      nextEdge == e1 ||
                      nextEdge == e2 ||
                      nextEdge == e3 ||
                      nextEdge == e4 )
               {
                  nextEdge++;
               }

               // Now we can safely collapse the edge.
               mesh.collapseEdge( e );
            }

            e = nextEdge;
         }

         cerr << "Improving vertex degrees..." << endl;

         // Next, we flip edges in an effort to get more uniform vertex valence.
         e = mesh.edgesBegin();
         while( e != mesh.edgesEnd() )
         {
            // After flipping, the current edge may no longer exist; therefore,
            // we need to grab a pointer to the next edge NOW.  (In principle,
            // one could implement edge flip such that the flipped edge is
            // not destroyed; however, we will not make that assumption here
            // for the sake of safety/robustness to future changes in the method
            // that performs the flipping.)
            EdgeIter nextEdge = e;
            nextEdge++;

            int a0 = e->halfedge()->vertex()->degree();
            int a1 = e->halfedge()->twin()->vertex()->degree();
            int b0 = e->halfedge()->next()->next()->vertex()->degree();
            int b1 = e->halfedge()->twin()->next()->next()->vertex()->degree();
            if( a0 != 3 && a1 != 3 && b0 != 3 && b1 != 3 )
            {
               int oldDefect = abs( a0-6 ) + abs( a1-6 ) + abs( b0-6 ) + abs( b1-6 );
               int newDefect = abs( (a0-1)-6 ) + abs( (a1-1)-6 ) + abs( (b0+1)-6 ) + abs( (b1+1)-6 );

               if( newDefect < oldDefect && a0-1 > 3 && a1-1 > 3 )
               {
                  mesh.flipEdge( e );
               }
            }

            e = nextEdge;
         }

         cerr << "Smoothing..." << endl;

         // Finally, we apply some tangential smoothing.
         for( int i = 0; i < nSmoothingIterations; i++ )
         {
            for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
            {
               v->computeCentroid();
            }

            for( VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ )
            {
               Vector3D N = v->normal();
               Vector3D c = v->centroid;
               Vector3D p = v->position;

               Vector3D u = 2*(c-p);
               u -= dot( u, N )*N;

               v->position += u;
            }
         }

         cerr << "Done!" << endl;
      }


   }
}
