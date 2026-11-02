/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2021-2024, INRIA.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** \author Jia Pan */

#ifndef COAL_TRAVERSAL_NODE_HFIELD_SHAPE_H
#define COAL_TRAVERSAL_NODE_HFIELD_SHAPE_H

#include <cmath>
#include <limits>

/// @cond INTERNAL

#include "coal/collision_data.h"
#include "coal/contact_patch.h"
#include "coal/shape/geometric_shapes.h"
#include "coal/narrowphase/narrowphase.h"
#include "coal/shape/geometric_shapes_utility.h"
#include "coal/internal/shape_shape_func.h"
#include "coal/internal/traversal_node_base.h"
#include "coal/internal/traversal.h"
#include "coal/internal/intersect.h"
#include "coal/hfield.h"
#include "coal/shape/convex.h"

namespace coal {

/// @addtogroup Traversal_For_Collision
/// @{

namespace details {
template <typename BV>
ConvexTpl<Quadrilateral32> buildConvexQuadrilateral(
    const HFNode<BV>& node, const HeightField<BV>& model) {
  const MatrixXs& heights = model.getHeights();
  const VecXs& x_grid = model.getXGrid();
  const VecXs& y_grid = model.getYGrid();

  const Scalar min_height = model.getMinHeight();

  const Scalar x0 = x_grid[node.x_id], x1 = x_grid[node.x_id + 1],
               y0 = y_grid[node.y_id], y1 = y_grid[node.y_id + 1];
  const Eigen::Block<const MatrixXs, 2, 2> cell =
      heights.block<2, 2>(node.y_id, node.x_id);

  assert(cell.maxCoeff() > min_height &&
         "max_height is lower than min_height");  // Check whether the geometry
                                                  // is degenerated

  std::shared_ptr<std::vector<Vec3s>> pts(new std::vector<Vec3s>({
      Vec3s(x0, y0, min_height),
      Vec3s(x0, y1, min_height),
      Vec3s(x1, y1, min_height),
      Vec3s(x1, y0, min_height),
      Vec3s(x0, y0, cell(0, 0)),
      Vec3s(x0, y1, cell(1, 0)),
      Vec3s(x1, y1, cell(1, 1)),
      Vec3s(x1, y0, cell(0, 1)),
  }));

  std::shared_ptr<std::vector<Quadrilateral32>> polygons(
      new std::vector<Quadrilateral32>(6));
  (*polygons)[0].set(0, 3, 2, 1);  // x+ side
  (*polygons)[1].set(0, 1, 5, 4);  // y- side
  (*polygons)[2].set(1, 2, 6, 5);  // x- side
  (*polygons)[3].set(2, 3, 7, 6);  // y+ side
  (*polygons)[4].set(3, 0, 4, 7);  // z- side
  (*polygons)[5].set(4, 5, 6, 7);  // z+ side

  return ConvexTpl<Quadrilateral32>(pts,  // points
                                    8,    // num points
                                    polygons,
                                    6  // number of polygons
  );
}

enum class FaceOrientationConvexPart1 {
  BOTTOM = 0,
  TOP = 1,
  WEST = 2,
  SOUTH_EAST = 4,
  NORTH = 8,
};

enum class FaceOrientationConvexPart2 {
  BOTTOM = 0,
  TOP = 1,
  SOUTH = 2,
  NORTH_WEST = 4,
  EAST = 8,
};

template <typename BV>
void buildConvexTriangles(const HFNode<BV>& node, const HeightField<BV>& model,
                          ConvexTpl<Triangle32>& convex1,
                          int& convex1_active_faces,
                          ConvexTpl<Triangle32>& convex2,
                          int& convex2_active_faces) {
  const MatrixXs& heights = model.getHeights();
  const VecXs& x_grid = model.getXGrid();
  const VecXs& y_grid = model.getYGrid();

  const Scalar min_height = model.getMinHeight();

  const Scalar x0 = x_grid[node.x_id], x1 = x_grid[node.x_id + 1],
               y0 = y_grid[node.y_id], y1 = y_grid[node.y_id + 1];
  const Scalar max_height = node.max_height;
  const Eigen::Block<const MatrixXs, 2, 2> cell =
      heights.block<2, 2>(node.y_id, node.x_id);

  const int contact_active_faces = node.contact_active_faces;
  convex1_active_faces = 0;
  convex2_active_faces = 0;

  typedef HFNodeBase::FaceOrientation FaceOrientation;

  if (contact_active_faces & FaceOrientation::TOP) {
    convex1_active_faces |= int(details::FaceOrientationConvexPart1::TOP);
    convex2_active_faces |= int(details::FaceOrientationConvexPart2::TOP);
  }

  if (contact_active_faces & FaceOrientation::BOTTOM) {
    convex1_active_faces |= int(details::FaceOrientationConvexPart1::BOTTOM);
    convex2_active_faces |= int(details::FaceOrientationConvexPart2::BOTTOM);
  }

  // Specific orientation for Convex1
  if (contact_active_faces & FaceOrientation::WEST) {
    convex1_active_faces |= int(details::FaceOrientationConvexPart1::WEST);
  }

  if (contact_active_faces & FaceOrientation::NORTH) {
    convex1_active_faces |= int(details::FaceOrientationConvexPart1::NORTH);
  }

  // Specific orientation for Convex2
  if (contact_active_faces & FaceOrientation::EAST) {
    convex2_active_faces |= int(details::FaceOrientationConvexPart2::EAST);
  }

  if (contact_active_faces & FaceOrientation::SOUTH) {
    convex2_active_faces |= int(details::FaceOrientationConvexPart2::SOUTH);
  }

  assert(max_height > min_height &&
         "max_height is lower than min_height");  // Check whether the geometry
                                                  // is degenerated
  COAL_UNUSED_VARIABLE(max_height);

  {
    std::shared_ptr<std::vector<Vec3s>> pts(new std::vector<Vec3s>({
        Vec3s(x0, y0, min_height),  // A
        Vec3s(x0, y1, min_height),  // B
        Vec3s(x1, y0, min_height),  // C
        Vec3s(x0, y0, cell(0, 0)),  // D
        Vec3s(x0, y1, cell(1, 0)),  // E
        Vec3s(x1, y0, cell(0, 1)),  // F
    }));

    std::shared_ptr<std::vector<Triangle32>> triangles(
        new std::vector<Triangle32>(8));
    (*triangles)[0].set(0, 2, 1);  // bottom
    (*triangles)[1].set(3, 4, 5);  // top
    (*triangles)[2].set(0, 1, 3);  // West 1
    (*triangles)[3].set(3, 1, 4);  // West 2
    (*triangles)[4].set(1, 2, 5);  // South-East 1
    (*triangles)[5].set(1, 5, 4);  // South-East 1
    (*triangles)[6].set(0, 5, 2);  // North 1
    (*triangles)[7].set(5, 0, 3);  // North 2

    convex1.set(pts,  // points
                6,    // num points
                triangles,
                8  // number of polygons
    );
  }

  {
    std::shared_ptr<std::vector<Vec3s>> pts(new std::vector<Vec3s>({
        Vec3s(x0, y1, min_height),  // A
        Vec3s(x1, y1, min_height),  // B
        Vec3s(x1, y0, min_height),  // C
        Vec3s(x0, y1, cell(1, 0)),  // D
        Vec3s(x1, y1, cell(1, 1)),  // E
        Vec3s(x1, y0, cell(0, 1)),  // F
    }));

    std::shared_ptr<std::vector<Triangle32>> triangles(
        new std::vector<Triangle32>(8));
    (*triangles)[0].set(2, 1, 0);  // bottom
    (*triangles)[1].set(3, 4, 5);  // top
    (*triangles)[2].set(0, 1, 3);  // South 1
    (*triangles)[3].set(3, 1, 4);  // South 2
    (*triangles)[4].set(0, 5, 2);  // North West 1
    (*triangles)[5].set(0, 3, 5);  // North West 2
    (*triangles)[6].set(1, 2, 5);  // East 1
    (*triangles)[7].set(4, 1, 2);  // East 2

    convex2.set(pts,  // points
                6,    // num points
                triangles,
                8  // number of polygons
    );
  }
}

inline Vec3s projectTriangle(const Vec3s& pointA, const Vec3s& pointB,
                             const Vec3s& pointC, const Vec3s& point) {
  const Project<Scalar>::ProjectResult result =
      Project<Scalar>::projectTriangle(pointA, pointB, pointC, point);
  Vec3s res = result.parameterization[0] * pointA +
              result.parameterization[1] * pointB +
              result.parameterization[2] * pointC;

  return res;
}

inline Vec3s projectTetrahedra(const Vec3s& pointA, const Vec3s& pointB,
                               const Vec3s& pointC, const Vec3s& pointD,
                               const Vec3s& point) {
  const Project<Scalar>::ProjectResult result =
      Project<Scalar>::projectTetrahedra(pointA, pointB, pointC, pointD, point);
  Vec3s res = result.parameterization[0] * pointA +
              result.parameterization[1] * pointB +
              result.parameterization[2] * pointC +
              result.parameterization[3] * pointD;

  return res;
}

inline Vec3s computeTriangleNormal(const Triangle32& triangle,
                                   const std::vector<Vec3s>& points) {
  const Vec3s pointA = points[triangle[0]];
  const Vec3s pointB = points[triangle[1]];
  const Vec3s pointC = points[triangle[2]];

  const Vec3s normal = (pointB - pointA).cross(pointC - pointA).normalized();
  assert(!normal.array().isNaN().any() && "normal is ill-defined");

  return normal;
}

inline Vec3s projectPointOnTriangle(const Vec3s& contact_point,
                                    const Triangle32& triangle,
                                    const std::vector<Vec3s>& points) {
  const Vec3s pointA = points[triangle[0]];
  const Vec3s pointB = points[triangle[1]];
  const Vec3s pointC = points[triangle[2]];

  const Vec3s contact_point_projected =
      projectTriangle(pointA, pointB, pointC, contact_point);

  return contact_point_projected;
}

inline Scalar distanceContactPointToTriangle(const Vec3s& contact_point,
                                             const Triangle32& triangle,
                                             const std::vector<Vec3s>& points) {
  const Vec3s contact_point_projected =
      projectPointOnTriangle(contact_point, triangle, points);
  return (contact_point_projected - contact_point).norm();
}

inline Scalar distanceContactPointToFace(const size_t face_id,
                                         const Vec3s& contact_point,
                                         const ConvexTpl<Triangle32>& convex,
                                         size_t& closest_face_id) {
  assert((face_id >= 0 && face_id < 8) && "face_id should be in [0;7]");

  const std::vector<Vec3s>& points = *(convex.points);
  if (face_id <= 1) {
    const Triangle32& triangle = (*(convex.polygons))[face_id];
    closest_face_id = face_id;
    return distanceContactPointToTriangle(contact_point, triangle, points);
  } else {
    const Triangle32& triangle1 = (*(convex.polygons))[face_id];
    const Scalar distance_to_triangle1 =
        distanceContactPointToTriangle(contact_point, triangle1, points);

    const Triangle32& triangle2 = (*(convex.polygons))[face_id + 1];
    const Scalar distance_to_triangle2 =
        distanceContactPointToTriangle(contact_point, triangle2, points);

    if (distance_to_triangle1 > distance_to_triangle2) {
      closest_face_id = face_id + 1;
      return distance_to_triangle2;
    } else {
      closest_face_id = face_id;
      return distance_to_triangle1;
    }
  }
}

template <typename Polygone, typename Shape>
bool binCorrection(const ConvexTpl<Polygone>& convex,
                   const int convex_active_faces, const Transform3s& tf1, const Shape& shape,
                   const Transform3s& shape_pose, Scalar& distance,
                   Vec3s& contact_1, Vec3s& contact_2, Vec3s& normal,
                   Vec3s& face_normal, const bool is_collision) {
  const Scalar prec = Scalar(1e-12);
  const std::vector<Vec3s>& points = *(convex.points);
  Vec3s contact_1_local = tf1.inverseTransform(contact_1);

  bool hfield_witness_is_on_bin_side = true;

  //  int closest_face_id_bottom_face = -1;
  //  int closest_face_id_top_face = -1;

  std::vector<size_t> active_faces;
  active_faces.reserve(5);
  active_faces.push_back(0);
  active_faces.push_back(1);

  if (convex_active_faces & 2) active_faces.push_back(2);
  if (convex_active_faces & 4) active_faces.push_back(4);
  if (convex_active_faces & 8) active_faces.push_back(6);

  Triangle32 face_triangle;
  Scalar shortest_distance_to_face = (std::numeric_limits<Scalar>::max)();
  face_normal = normal;
  for (const size_t active_face : active_faces) {
    size_t closest_face_id;
    const Scalar distance_to_face = distanceContactPointToFace(
        active_face, contact_1_local, convex, closest_face_id);

    const bool contact_point_is_on_face = distance_to_face <= prec;
    if (contact_point_is_on_face) {
      hfield_witness_is_on_bin_side = false;
      face_triangle = (*(convex.polygons))[closest_face_id];
      shortest_distance_to_face = distance_to_face;
      break;
    } else if (distance_to_face < shortest_distance_to_face) {
      face_triangle = (*(convex.polygons))[closest_face_id];
      shortest_distance_to_face = distance_to_face;
    }
  }

  // We correct only if there is a collision with the bin
  if (is_collision) {
    if (!face_triangle.isValid())
      COAL_THROW_PRETTY("face_triangle is not initialized", std::logic_error);

    const Vec3s face_pointA = points[face_triangle[0]];
    // face_normal = computeTriangleNormal(face_triangle, points);
    Vec3s face_normal_local = computeTriangleNormal(face_triangle, points);
    Transform3s shape_pose_local = tf1.inverseTimes(shape_pose);
    int hint = 0;
    // Since we compute the support manually, we need to take into account the
    // sphere swept radius of the shape.
    // TODO: take into account the swept-sphere radius of the bin.
    const Vec3s _support = getSupport<details::SupportOptions::WithSweptSphere>(
        &shape, -shape_pose_local.rotation().transpose() * face_normal_local, hint);
    const Vec3s support_local =
        shape_pose_local.rotation() * _support + shape_pose_local.translation();

    // Project support into the inclined bin having triangle
    const Scalar offset_plane = face_normal_local.dot(face_pointA);
    const Plane projection_plane(face_normal_local, offset_plane);
    const Scalar distance_support_projection_plane = projection_plane.signedDistance(support_local);

    Vec3s projected_support_local = support_local - distance_support_projection_plane * face_normal_local;
    Vec3s contact_1_out_local = projectPointOnTriangle(projected_support_local, face_triangle, points);
    Vec3s contact_2_out_local = contact_1_out_local + distance_support_projection_plane * face_normal_local;

    contact_1 = tf1.transform(contact_1_out_local);
    contact_2 = tf1.transform(contact_2_out_local);
    normal = tf1.getRotation() * face_normal_local;
    face_normal = normal;
    distance = -std::fabs(distance_support_projection_plane);
  }

  return hfield_witness_is_on_bin_side;
}

template <typename Polygone, typename Shape, int Options>
bool shapeDistance(const GJKSolver* nsolver, const CollisionRequest& request,
                   const ConvexTpl<Polygone>& convex1,
                   const int convex1_active_faces,
                   const ConvexTpl<Polygone>& convex2,
                   const int convex2_active_faces, const Transform3s& tf1,
                   const Shape& shape, const Transform3s& tf2, Scalar& distance,
                   Vec3s& c1, Vec3s& c2, Vec3s& normal, Vec3s& normal_top,
                   bool& hfield_witness_is_on_bin_side) {
  enum { RTIsIdentity = Options & RelativeTransformationIsIdentity };

  const Transform3s Id;
  // The solver `nsolver` has already been set up by the collision request
  // `request`. If GJK early stopping is enabled through `request`, it will be
  // used.
  // The only thing we need to make sure is that in case of collision, the
  // penetration information is computed (as we do bins comparison).
  const bool compute_penetration = true;
  Vec3s contact1_1, contact1_2, contact2_1, contact2_2;
  Vec3s normal1, normal1_top, normal2, normal2_top;
  Scalar distance1, distance2;

  if (RTIsIdentity) {
    distance1 = internal::ShapeShapeDistance<ConvexTpl<Polygone>, Shape>(
        &convex1, Id, &shape, tf2, nsolver, compute_penetration, contact1_1,
        contact1_2, normal1);
  } else {
    distance1 = internal::ShapeShapeDistance<ConvexTpl<Polygone>, Shape>(
        &convex1, tf1, &shape, tf2, nsolver, compute_penetration, contact1_1,
        contact1_2, normal1);
  }
  bool collision1 = (distance1 - request.security_margin <=
                     request.collision_distance_threshold);

  bool hfield_witness_is_on_bin_side1 =
      binCorrection(convex1, convex1_active_faces, tf1, shape, tf2, distance1,
                    contact1_1, contact1_2, normal1, normal1_top, collision1);

  if (RTIsIdentity) {
    distance2 = internal::ShapeShapeDistance<ConvexTpl<Polygone>, Shape>(
        &convex2, Id, &shape, tf2, nsolver, compute_penetration, contact2_1,
        contact2_2, normal2);
  } else {
    distance2 = internal::ShapeShapeDistance<ConvexTpl<Polygone>, Shape>(
        &convex2, tf1, &shape, tf2, nsolver, compute_penetration, contact2_1,
        contact2_2, normal2);
  }
  bool collision2 = (distance2 - request.security_margin <=
                     request.collision_distance_threshold);

  bool hfield_witness_is_on_bin_side2 =
      binCorrection(convex2, convex2_active_faces, tf1, shape, tf2, distance2,
                    contact2_1, contact2_2, normal2, normal2_top, collision2);

  if (collision1 && collision2) {
    if (distance1 > distance2)  // switch values
    {
      distance = distance2;
      c1 = contact2_1;
      c2 = contact2_2;
      normal = normal2;
      normal_top = normal2_top;
      hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
    } else {
      distance = distance1;
      c1 = contact1_1;
      c2 = contact1_2;
      normal = normal1;
      normal_top = normal1_top;
      hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
    }
    return true;
  } else if (collision1) {
    distance = distance1;
    c1 = contact1_1;
    c2 = contact1_2;
    normal = normal1;
    normal_top = normal1_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
    return true;
  } else if (collision2) {
    distance = distance2;
    c1 = contact2_1;
    c2 = contact2_2;
    normal = normal2;
    normal_top = normal2_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
    return true;
  }

  if (distance1 > distance2)  // switch values
  {
    distance = distance2;
    c1 = contact2_1;
    c2 = contact2_2;
    normal = normal2;
    normal_top = normal2_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
  } else {
    distance = distance1;
    c1 = contact1_1;
    c2 = contact1_2;
    normal = normal1;
    normal_top = normal1_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
  }
  return false;
}

}  // namespace details

/// @brief Traversal node for collision between height field and shape
template <typename BV, typename S,
          int _Options = RelativeTransformationIsIdentity>
class HeightFieldShapeCollisionTraversalNode
    : public CollisionTraversalNodeBase {
 public:
  typedef CollisionTraversalNodeBase Base;

  enum {
    Options = _Options,
    RTIsIdentity = _Options & RelativeTransformationIsIdentity
  };

  HeightFieldShapeCollisionTraversalNode(const CollisionRequest& request)
      : CollisionTraversalNodeBase(request) {
    model1 = NULL;
    model2 = NULL;

    num_bv_tests = 0;
    num_leaf_tests = 0;
    query_time_seconds = 0.0;

    nsolver = NULL;
    count = 0;
  }

  /// @brief Whether the BV node in the first BVH tree is leaf
  bool isFirstNodeLeaf(unsigned int b) const {
    return model1->getBV(b).isLeaf();
  }

  /// @brief Obtain the left child of BV node in the first BVH
  int getFirstLeftChild(unsigned int b) const {
    return static_cast<int>(model1->getBV(b).leftChild());
  }

  /// @brief Obtain the right child of BV node in the first BVH
  int getFirstRightChild(unsigned int b) const {
    return static_cast<int>(model1->getBV(b).rightChild());
  }

  /// test between BV b1 and shape
  /// @param b1 BV to test,
  /// @retval sqrDistLowerBound square of a lower bound of the minimal
  ///         distance between bounding volumes.
  /// @brief BV culling test in one BVTT node
  bool BVDisjoints(unsigned int b1, unsigned int /*b2*/,
                   Scalar& sqrDistLowerBound) const {
    if (this->enable_statistics) this->num_bv_tests++;

    bool disjoint;
    if (RTIsIdentity) {
      assert(false && "must never happened");
      disjoint = !this->model1->getBV(b1).bv.overlap(
          this->model2_bv, this->request, sqrDistLowerBound);
    } else {
      disjoint = !overlap(this->tf1.getRotation(), this->tf1.getTranslation(),
                          this->model1->getBV(b1).bv, this->model2_bv,
                          this->request, sqrDistLowerBound);
    }

    if (disjoint)
      internal::updateDistanceLowerBoundFromBV(this->request, *this->result,
                                               sqrDistLowerBound);

    assert(!disjoint || sqrDistLowerBound > 0);
    return disjoint;
  }

  /// @brief Intersection testing between leaves (one ConvexTpl and one shape)
  void leafCollides(unsigned int b1, unsigned int /*b2*/,
                    Scalar& sqrDistLowerBound) const {
    count++;
    if (this->enable_statistics) this->num_leaf_tests++;
    const HFNode<BV>& node = this->model1->getBV(b1);

    // Split quadrilateral primitives into two convex shapes corresponding to
    // polyhedron with triangular bases. This is essential to keep the convexity

    //    typedef ConvexTpl<Quadrilateral32> ConvexQuadrilateral32;
    //    const ConvexQuadrilateral32 convex =
    //    details::buildConvexQuadrilateral(node,*this->model1);

    typedef ConvexTpl<Triangle32> ConvexTriangle32;
    ConvexTriangle32 convex1, convex2;
    int convex1_active_faces, convex2_active_faces;
    // TODO: inherit from hfield's inflation here
    details::buildConvexTriangles(node, *this->model1, convex1,
                                  convex1_active_faces, convex2,
                                  convex2_active_faces);

    // Compute aabb_local for BoundingVolumeGuess case in the GJK solver
    if (nsolver->gjk_initial_guess == GJKInitialGuess::BoundingVolumeGuess) {
      convex1.computeLocalAABB();
      convex2.computeLocalAABB();
    }

    Scalar distance;
    //    Vec3s contact_point, normal;
    Vec3s c1, c2, normal, normal_face;
    bool hfield_witness_is_on_bin_side;

    bool collision = details::shapeDistance<Triangle32, S, Options>(
        nsolver, this->request, convex1, convex1_active_faces, convex2,
        convex2_active_faces, this->tf1, *(this->model2), this->tf2, distance,
        c1, c2, normal, normal_face, hfield_witness_is_on_bin_side);

    Scalar distToCollision = distance - this->request.security_margin;
    if (distToCollision <= this->request.collision_distance_threshold) {
      sqrDistLowerBound = 0;
      if (this->result->numContacts() < this->request.num_max_contacts) {
        if (normal_face.isApprox(normal) &&
            (collision || !hfield_witness_is_on_bin_side)) {
          const Contact leaf_contact(this->model1, this->model2, (int)b1,
                                     (int)Contact::NONE, c1, c2, normal,
                                     distance);
          const NODE_TYPE shape_type = this->model2->getNodeType();
          const auto cylinder_witness_is_valid = [&](const Vec3s& witness) {
            const Cylinder* cylinder =
                dynamic_cast<const Cylinder*>(this->model2);
            assert(cylinder != nullptr);
            const Vec3s local = this->tf2.inverseTransform(witness);
            const Scalar radial = local.head<2>().norm();
            const Scalar axial = std::abs(local.z());
            const Scalar tolerance = Scalar(1e-5);
            return radial <= cylinder->radius + tolerance &&
                   axial <= cylinder->halfLength + tolerance &&
                   (std::abs(radial - cylinder->radius) <= tolerance ||
                    std::abs(axial - cylinder->halfLength) <= tolerance);
          };
          const auto cylinder_surface_witness =
              [&](const Vec3s& terrain_witness, const Vec3s& contact_normal,
                  Vec3s& shape_witness, Scalar& contact_depth) {
                const Cylinder* cylinder =
                    dynamic_cast<const Cylinder*>(this->model2);
                assert(cylinder != nullptr);
                const Vec3s origin =
                    this->tf2.inverseTransform(terrain_witness);
                const Vec3s direction =
                    this->tf2.rotation().transpose() * contact_normal;
                const Scalar radius2 = cylinder->radius * cylinder->radius;
                Scalar best_t = -std::numeric_limits<Scalar>::infinity();
                const auto accept = [&](Scalar t) {
                  if (t <= Scalar(1e-6) && t > best_t) best_t = t;
                };
                if (std::abs(direction.z()) > Scalar(1e-8)) {
                  for (Scalar cap_z : {-cylinder->halfLength,
                                       cylinder->halfLength}) {
                    const Scalar t = (cap_z - origin.z()) / direction.z();
                    const Vec3s point = origin + t * direction;
                    if (point.head<2>().squaredNorm() <=
                        radius2 + Scalar(1e-8))
                      accept(t);
                  }
                }
                const Scalar a = direction.head<2>().squaredNorm();
                if (a > Scalar(1e-12)) {
                  const Scalar b = Scalar(2) *
                      origin.head<2>().dot(direction.head<2>());
                  const Scalar c = origin.head<2>().squaredNorm() - radius2;
                  const Scalar discriminant = b * b - Scalar(4) * a * c;
                  if (discriminant >= Scalar(0)) {
                    const Scalar root = std::sqrt(discriminant);
                    for (Scalar t : {(-b - root) / (Scalar(2) * a),
                                     (-b + root) / (Scalar(2) * a)}) {
                      const Scalar z = origin.z() + t * direction.z();
                      if (std::abs(z) <=
                          cylinder->halfLength + Scalar(1e-8))
                        accept(t);
                    }
                  }
                }
                if (!std::isfinite(best_t)) return false;
                shape_witness =
                    this->tf2.transform(origin + best_t * direction);
                contact_depth = best_t;
                return true;
              };
          bool generated_patch_contact = false;
          if (this->request.enable_height_field_contact_patches) {
            const auto& cell = this->model1->getBV(b1);
            const VecXs& x_grid = this->model1->getXGrid();
            const VecXs& y_grid = this->model1->getYGrid();
            const MatrixXs& heights = this->model1->getHeights();
            const Scalar x0 = x_grid[cell.x_id];
            const Scalar x1 = x_grid[cell.x_id + 1];
            const Scalar y0 = y_grid[cell.y_id];
            const Scalar y1 = y_grid[cell.y_id + 1];
            const auto surface_height = [&](const Scalar u, const Scalar v) {
              const Scalar h00 = heights(cell.y_id, cell.x_id);
              const Scalar h01 = heights(cell.y_id + 1, cell.x_id);
              const Scalar h10 = heights(cell.y_id, cell.x_id + 1);
              const Scalar h11 = heights(cell.y_id + 1, cell.x_id + 1);
              return u + v <= Scalar(1)
                         ? h00 + u * (h10 - h00) + v * (h01 - h00)
                         : h11 + (Scalar(1) - u) * (h01 - h11) +
                               (Scalar(1) - v) * (h10 - h11);
            };
            if (shape_type == GEOM_SPHERE || shape_type == GEOM_CAPSULE) {
              // Curved shape patches must use points on the actual surface.
              // Leaf witnesses from neighbouring cells are not substitutes
              // for the support points of this face.
              generated_patch_contact = true;
              const auto curved_supports = [&](const Vec3s& world_normal) {
                const Vec3s shape_normal =
                    this->tf2.rotation().transpose() * world_normal;
                std::vector<Vec3s> supports;
                if (shape_type == GEOM_SPHERE) {
                  const Sphere* sphere =
                      dynamic_cast<const Sphere*>(this->model2);
                  assert(sphere != nullptr);
                  supports.push_back(-sphere->radius * shape_normal);
                } else {
                  const Capsule* capsule =
                      dynamic_cast<const Capsule*>(this->model2);
                  assert(capsule != nullptr);
                  Vec3s support = -capsule->radius * shape_normal;
                  if (std::abs(shape_normal.z()) < Scalar(1e-3)) {
                    support.z() = -capsule->halfLength;
                    supports.push_back(support);
                    support.z() = capsule->halfLength;
                    supports.push_back(support);
                  } else {
                    support.z() += shape_normal.z() > 0
                                       ? -capsule->halfLength
                                       : capsule->halfLength;
                    supports.push_back(support);
                  }
                }
                return supports;
              };
              for (const Vec3s& shape_local : curved_supports(normal)) {
                const Vec3s shape_witness = this->tf2.transform(shape_local);
                const Vec3s shape_terrain_local =
                    this->tf1.inverseTransform(shape_witness);
                const Scalar u0 = (shape_terrain_local.x() - x0) / (x1 - x0);
                const Scalar v0 = (y0 - shape_terrain_local.y()) / (y0 - y1);
                const Vec3s normal_terrain =
                    this->tf1.rotation().transpose() * normal;
                if (std::abs(normal_terrain.z()) < Scalar(1e-4)) continue;
                const Scalar h00 = heights(cell.y_id, cell.x_id);
                const Scalar h01 = heights(cell.y_id + 1, cell.x_id);
                const Scalar h10 = heights(cell.y_id, cell.x_id + 1);
                const Scalar h11 = heights(cell.y_id + 1, cell.x_id + 1);
                Scalar u_face = u0, v_face = v0;
                Vec3s terrain_local_surface;
                for (int iteration = 0; iteration < 2; ++iteration) {
                  const bool lower_triangle = u_face + v_face <= Scalar(1);
                  const Scalar dh_dx =
                      (lower_triangle ? h10 - h00 : h11 - h01) / (x1 - x0);
                  const Scalar dh_dy =
                      (lower_triangle ? h01 - h00 : h11 - h10) / (y1 - y0);
                  const Scalar denominator = normal_terrain.z() -
                      dh_dx * normal_terrain.x() -
                      dh_dy * normal_terrain.y();
                  if (std::abs(denominator) < Scalar(1e-6)) break;
                  const Scalar separation =
                      (shape_terrain_local.z() -
                       (lower_triangle
                            ? h00 + u0 * (h10 - h00) + v0 * (h01 - h00)
                            : h11 + (Scalar(1) - u0) * (h01 - h11) +
                                  (Scalar(1) - v0) * (h10 - h11))) /
                      denominator;
                  terrain_local_surface =
                      shape_terrain_local - separation * normal_terrain;
                  u_face = (terrain_local_surface.x() - x0) / (x1 - x0);
                  v_face = (y0 - terrain_local_surface.y()) / (y0 - y1);
                  if ((u_face + v_face <= Scalar(1)) == lower_triangle) break;
                }
                const Scalar u = (terrain_local_surface.x() - x0) / (x1 - x0);
                const Scalar v = (y0 - terrain_local_surface.y()) / (y0 - y1);
                if (u < Scalar(-1e-6) || u > Scalar(1) + Scalar(1e-6) ||
                    v < Scalar(-1e-6) || v > Scalar(1) + Scalar(1e-6) ||
                    std::abs(terrain_local_surface.z() - surface_height(u, v)) >
                        Scalar(1e-5))
                  continue;
                const Vec3s terrain_witness =
                    this->tf1.transform(terrain_local_surface);
                const Scalar depth =
                    (shape_witness - terrain_witness).dot(normal);
                const Contact curved_contact(
                    this->model1, this->model2, (int)b1,
                    (int)Contact::NONE, terrain_witness, shape_witness, normal,
                    depth);
                this->result->addContactIfDistinct(curved_contact);
              }
              const Vec3s center_terrain =
                  this->tf1.inverseTransform(this->tf2.translation());
              const auto append_wall = [&](int axis, Scalar wall,
                                           const Vec3s& wall_normal_local) {
                const Vec3s wall_normal =
                    this->tf1.rotation() * wall_normal_local;
                for (const Vec3s& shape_local : curved_supports(wall_normal)) {
                  const Vec3s shape_witness =
                      this->tf2.transform(shape_local);
                  const Vec3s shape_terrain =
                      this->tf1.inverseTransform(shape_witness);
                  Vec3s wall_point = shape_terrain;
                  wall_point[axis] = wall;
                  const Scalar u = (wall_point.x() - x0) / (x1 - x0);
                  const Scalar v = (y0 - wall_point.y()) / (y0 - y1);
                  if (u < Scalar(-1e-6) || u > Scalar(1) + Scalar(1e-6) ||
                      v < Scalar(-1e-6) || v > Scalar(1) + Scalar(1e-6) ||
                      wall_point.z() < this->model1->getMinHeight() ||
                      wall_point.z() > surface_height(u, v))
                    continue;
                  const Vec3s terrain_witness = this->tf1.transform(wall_point);
                  const Scalar depth =
                      (shape_witness - terrain_witness).dot(wall_normal);
                  if (depth >= Scalar(-1e-8)) continue;
                  const Contact wall_contact(
                      this->model1, this->model2, (int)b1,
                      (int)Contact::NONE, terrain_witness, shape_witness,
                      wall_normal, depth);
                  this->result->addContactIfDistinct(wall_contact);
                }
              };
              if (cell.x_id == 0 && center_terrain.x() < x_grid[0])
                append_wall(0, x_grid[0], Vec3s(-1, 0, 0));
              if (cell.x_id + 2 == x_grid.size() &&
                  center_terrain.x() > x_grid[x_grid.size() - 1])
                append_wall(0, x_grid[x_grid.size() - 1], Vec3s(1, 0, 0));
              if (cell.y_id == 0 && center_terrain.y() > y_grid[0])
                append_wall(1, y_grid[0], Vec3s(0, 1, 0));
              if (cell.y_id + 2 == y_grid.size() &&
                  center_terrain.y() < y_grid[y_grid.size() - 1])
                append_wall(1, y_grid[y_grid.size() - 1], Vec3s(0, -1, 0));
            } else {
              CollisionResult witness_result;
              witness_result.addContact(leaf_contact);
              const ContactPatchRequest patch_request(1, 16, Scalar(1e-3));
              const auto append_patch = [&](const CollisionGeometry* prism) {
              ContactPatchResult patch_result(patch_request);
              computeContactPatch(prism, this->tf1, this->model2, this->tf2,
                                  witness_result, patch_request,
                                  patch_result);
              for (size_t patch_index = 0;
                   patch_index < patch_result.numContactPatches();
                   ++patch_index) {
                const ContactPatch& patch =
                    patch_result.getContactPatch(patch_index);
                const Vec3s patch_normal = patch.getNormal();
                const Scalar patch_depth = patch.penetration_depth;
                for (size_t point_index = 0; point_index < patch.size();
                     ++point_index) {
                  const Vec3s midpoint = patch.getPoint(point_index);
                  const Vec3s terrain_witness =
                      midpoint - Scalar(0.5) * patch_depth * patch_normal;
                  const Vec3s terrain_local =
                      this->tf1.inverseTransform(terrain_witness);
                  const Scalar u = (terrain_local.x() - x0) / (x1 - x0);
                  const Scalar v = (y0 - terrain_local.y()) / (y0 - y1);
                  if (u < Scalar(-1e-6) || u > Scalar(1) + Scalar(1e-6) ||
                      v < Scalar(-1e-6) || v > Scalar(1) + Scalar(1e-6))
                    continue;
                  if (std::abs(terrain_local.z() - surface_height(u, v)) >
                      Scalar(1e-5))
                    continue;
                  Vec3s shape_witness =
                      midpoint + Scalar(0.5) * patch_depth * patch_normal;
                  Scalar contact_depth = patch_depth;
                  if (shape_type == GEOM_CYLINDER &&
                      !cylinder_witness_is_valid(shape_witness)) {
                    if (!cylinder_surface_witness(
                            terrain_witness, patch_normal, shape_witness,
                            contact_depth))
                      continue;
                  }
                  const Contact patch_contact(
                      this->model1, this->model2, (int)b1,
                      (int)Contact::NONE, terrain_witness, shape_witness,
                      patch_normal, contact_depth);
                  this->result->addContactIfDistinct(patch_contact);
                  generated_patch_contact = true;
                }
              }
              };
              append_patch(&convex1);
              append_patch(&convex2);
            }
          }
          if (!generated_patch_contact) {
            if (this->request.enable_height_field_contact_patches &&
                shape_type == GEOM_CYLINDER &&
                !cylinder_witness_is_valid(c2)) {
              Vec3s shape_witness;
              Scalar contact_depth;
              if (cylinder_surface_witness(c1, normal, shape_witness,
                                           contact_depth)) {
                const Contact cylinder_contact(
                    this->model1, this->model2, (int)b1,
                    (int)Contact::NONE, c1, shape_witness, normal,
                    contact_depth);
                this->result->addContactIfDistinct(cylinder_contact);
              }
            } else {
              this->result->addContact(leaf_contact);
            }
          }
          assert(this->result->isCollision());
        }
      }
    } else
      sqrDistLowerBound = distToCollision * distToCollision;

    //    const Vec3s c1 = contact_point - distance * 0.5 * normal;
    //    const Vec3s c2 = contact_point + distance * 0.5 * normal;
    internal::updateDistanceLowerBoundFromLeaf(this->request, *this->result,
                                               distToCollision, c1, c2, normal);

    assert(this->result->isCollision() || sqrDistLowerBound > 0);
  }

  const GJKSolver* nsolver;

  const HeightField<BV>* model1;
  const S* model2;
  BV model2_bv;

  mutable int num_bv_tests;
  mutable int num_leaf_tests;
  mutable Scalar query_time_seconds;
  mutable int count;
};

/// @}

/// @addtogroup Traversal_For_Distance
/// @{

/// @brief Traversal node for distance between height field and shape
template <typename BV, typename S,
          int _Options = RelativeTransformationIsIdentity>
class HeightFieldShapeDistanceTraversalNode : public DistanceTraversalNodeBase {
 public:
  typedef DistanceTraversalNodeBase Base;

  enum {
    Options = _Options,
    RTIsIdentity = _Options & RelativeTransformationIsIdentity
  };

  HeightFieldShapeDistanceTraversalNode() : DistanceTraversalNodeBase() {
    model1 = NULL;
    model2 = NULL;

    num_leaf_tests = 0;
    query_time_seconds = 0.0;

    rel_err = 0;
    abs_err = 0;
    nsolver = NULL;
  }

  /// @brief Whether the BV node in the first BVH tree is leaf
  bool isFirstNodeLeaf(unsigned int b) const {
    return model1->getBV(b).isLeaf();
  }

  /// @brief Obtain the left child of BV node in the first BVH
  int getFirstLeftChild(unsigned int b) const {
    return model1->getBV(b).leftChild();
  }

  /// @brief Obtain the right child of BV node in the first BVH
  int getFirstRightChild(unsigned int b) const {
    return model1->getBV(b).rightChild();
  }

  /// @brief BV culling test in one BVTT node
  Scalar BVDistanceLowerBound(unsigned int b1, unsigned int /*b2*/) const {
    return model1->getBV(b1).bv.distance(
        model2_bv);  // TODO(jcarpent): tf1 is not taken into account here.
  }

  /// @brief Distance testing between leaves (one bin of the height field and
  /// one shape)
  /// TODO(louis): deal with Hfield-Shape distance just like in Hfield-Shape
  /// collision (bin correction etc).
  void leafComputeDistance(unsigned int b1, unsigned int /*b2*/) const {
    if (this->enable_statistics) this->num_leaf_tests++;

    const BVNode<BV>& node = this->model1->getBV(b1);

    typedef ConvexTpl<Quadrilateral32> ConvexQuadrilateral32;
    const ConvexQuadrilateral32 convex =
        details::buildConvexQuadrilateral(node, *this->model1);

    Vec3s p1, p2, normal;
    const Scalar distance =
        internal::ShapeShapeDistance<ConvexQuadrilateral32, S>(
            &convex, this->tf1, this->model2, this->tf2, this->nsolver,
            this->request.enable_signed_distance, p1, p2, normal);

    this->result->update(distance, this->model1, this->model2, b1,
                         DistanceResult::NONE, p1, p2, normal);
  }

  /// @brief Whether the traversal process can stop early
  bool canStop(Scalar c) const {
    if ((c >= this->result->min_distance - abs_err) &&
        (c * (1 + rel_err) >= this->result->min_distance))
      return true;
    return false;
  }

  Scalar rel_err;
  Scalar abs_err;

  const GJKSolver* nsolver;

  const HeightField<BV>* model1;
  const S* model2;
  BV model2_bv;

  mutable int num_bv_tests;
  mutable int num_leaf_tests;
  mutable Scalar query_time_seconds;
};

/// @}

}  // namespace coal

/// @endcond

#endif
