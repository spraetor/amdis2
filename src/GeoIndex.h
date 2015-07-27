#pragma once

namespace AMDiS 
{
  /// internal used indices to represent the different geometrical objects.
  /// Used as parameter for getGeo() and as template parameter for FixVec. 
  enum GeoIndex
  {
    CENTER   = 0, /**< in 1d the center is at the edge, in 2d at the face, in 3d 
		    * at the interior of an element. So a Line has no edge but
		    * only a center, a Triangle has no face but only a center.
		    */
    VERTEX   = 1, /**< index for element vertices.
		    * number of vertices is equal to number of parts and 
		    * neighbours.
		    */
    EDGE     = 2, /**< index for element edges */
    FACE     = 3, /**< index for element faces */
    DIMEN    =-1, /**< index for problem dimension */
    PARTS    =-2, /**< index for parts of an element (vertices in 1d, edges in 2d
		    * , faces in 3d). Number of element parts is equal to number
		    * of vertices and neighbours. 
		    */
    NEIGH    =-3, /**< index for neighbours of an element.
		    * Number of element neighbours is equal to number of 
		    * vertices and parts.
		    */
    WORLD    =-4, /**< index for world dimension */
    BOUNDARY =-5, /**< index for boundary nodes of an element. This could be
		    * vertices, edges or faces.
		    */
    PROJECTION=-6, /**< index for element and boundary projections */
    
    NO_INDEX =-127
  };

  /// \cond HIDDEN_SYMBOLS
  // dummy-type for partial template specialization
  template <GeoIndex> struct _geo {};
  /// \endcond

#define MAXPART FACE
#define MINPART PROJECTION


  /// Returns the GeoIndex of d for dimension dim.
#define INDEX_OF_DIM(d, dim) (static_cast<GeoIndex>((d == dim) ? CENTER : d + 1))

  /// Returns the dimension of GeoIndex ind for dimension dim
#define DIM_OF_INDEX(ind, dim) ((static_cast<int>(ind) == 0) ? dim : static_cast<int>(ind) - 1)

} // end namespace AMDiS
